@testset "certificates on hand-made problems" begin
    lp = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [3, 2],
        [([1, 1], :le, 4), ([1, 3], :le, 6), ([1, 0], :le, 3)])
    s = solve(lp)
    @test verify(lp, s)
    # shadow prices: relaxing rows 1 and 3 helps, row 2 is slack
    @test s.duals == [2, 0, 1]
    @test s.reduced_costs == [0, 0]
    @test check_certificate(lp, s) === nothing

    # tampering with any part of the certificate is caught
    @test !verify(lp, Solution(:optimal; objective = s.objective, x = s.x, duals = [2, 0, 0], reduced_costs = s.reduced_costs))
    @test !verify(lp, Solution(:optimal; objective = s.objective - 1, x = s.x, duals = s.duals))
    @test !verify(lp, Solution(:optimal; objective = s.objective, x = Q[4, 0], duals = s.duals))
    @test !verify(lp, Solution(:optimal; objective = s.objective, x = s.x, duals = [-2, 0, 1]))
    suboptimal = Q[2, 1]
    @test !verify(lp, Solution(:optimal; objective = objective_value(lp, suboptimal), x = suboptimal, duals = s.duals))

    inf = lp_of(:min, [("x", 0, nothing), ("y", 0, nothing)], [1, 1],
        [([1, 1], :le, 1), ([1, 1], :ge, 3)])
    si = solve(inf)
    @test si.status === :infeasible
    @test verify(inf, si)
    @test !verify(inf, Solution(:infeasible; farkas = Q[0, 0]))
    @test !verify(inf, Solution(:infeasible; farkas = Q[1, 0]))
    @test !verify(inf, Solution(:infeasible; farkas = Q[-1, 1]))   # wrong signs

    unb = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [1, 1], [([1, -1], :le, 1)])
    su = solve(unb)
    @test su.status === :unbounded && verify(unb, su)
    @test !verify(unb, Solution(:unbounded; x = su.x, ray = Q[1, 1] .* 0))
    @test !verify(unb, Solution(:unbounded; x = su.x, ray = Q[1, 0]))     # leaves the feasible region
    @test !verify(unb, Solution(:unbounded; x = Q[-1, 0], ray = su.ray))  # infeasible start

    clash = lp_of(:min, [("x", 3, 2)], [1], [])
    @test verify(clash, solve(clash))
    @test verify(LP(), solve(LP()))
end

@testset "duals are shadow prices" begin
    # min 2x + 3y s.t. x + y >= 10 (row 1), x <= 8 (row 2) -> x = 8, y = 2
    lp = lp_of(:min, [("x", 0, nothing), ("y", 0, nothing)], [2, 3],
        [([1, 1], :ge, 10), ([1, 0], :le, 8)])
    s = solve(lp)
    @test s.objective == 22 && s.duals == [3, -1]
    for (i, delta) in ((1, 1), (2, 1))
        lp2 = lp_of(:min, [("x", 0, nothing), ("y", 0, nothing)], [2, 3],
            [([1, 1], :ge, 10 + (i == 1) * delta), ([1, 0], :le, 8 + (i == 2) * delta)])
        @test solve(lp2).objective - s.objective == s.duals[i] * delta
    end
    # a max problem reports d objective / d rhs as well (row 2 is slack here)
    mk(rhs1) = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [3, 2],
        [([1, 1], :le, rhs1), ([1, 3], :le, 9), ([1, 0], :le, 3)])
    sm = solve(mk(4))
    @test sm.duals == [2, 0, 1]
    @test solve(mk(5)).objective - sm.objective == sm.duals[1]
end

@testset "random LPs agree with brute-force vertex enumeration" begin
    rng = MersenneTwister(20261105)
    counts = Dict(:optimal => 0, :infeasible => 0)
    for t in 1:400
        lp = random_lp(rng; boxed = true)
        s = solve(lp)
        status, best = brute_force(lp)
        @test s.status === status
        if status === :optimal
            @test s.objective == best
        end
        @test check_certificate(lp, s) === nothing
        counts[s.status] += 1
    end
    @test counts[:optimal] > 100 && counts[:infeasible] > 20
end

@testset "random LPs with open bounds always produce valid certificates" begin
    rng = MersenneTwister(77)
    seen = Dict(:optimal => 0, :infeasible => 0, :unbounded => 0)
    for t in 1:800
        lp = random_lp(rng; boxed = false, maxn = 4, maxm = 4)
        s = solve(lp)
        reason = check_certificate(lp, s)
        @test reason === nothing
        reason === nothing || println(reason, " in ", lp)
        for rule in (:bland, :dantzig)
            t > 300 && break
            s2 = solve(lp; rule)
            @test s2.status === s.status
            s.status === :optimal && @test s2.objective == s.objective
            @test verify(lp, s2)
        end
        seen[s.status] += 1
    end
    @test all(>(30), values(seen))
end
