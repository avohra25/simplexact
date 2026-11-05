function knapsack(values, weights, cap)
    lp = LP(:max)
    xs = [add_variable!(lp, "x$i"; ub = 1, integer = true) for i in eachindex(values)]
    set_objective!(lp, [xs[i] => values[i] for i in eachindex(xs)])
    add_constraint!(lp, [xs[i] => weights[i] for i in eachindex(xs)], :le, cap)
    lp, xs
end

function knapsack_brute(values, weights, cap)
    best = 0
    n = length(values)
    for mask in 0:(2^n-1)
        w = sum((weights[i] for i in 1:n if (mask >> (i - 1)) & 1 == 1); init = 0)
        w <= cap || continue
        best = max(best, sum((values[i] for i in 1:n if (mask >> (i - 1)) & 1 == 1); init = 0))
    end
    best
end

"Exhaustive oracle for MIPs with boxed integer variables: try every integer
assignment and solve the remaining LP over the continuous variables."
function mip_brute(lp::LP)
    ints = [j for j in 1:nvars(lp) if lp.integer[j]]
    lb, ub = Simplexact.tightened_bounds(lp)
    ranges = [Int(lb[j]):Int(ub[j]) for j in ints]
    best = nothing
    s = lp.sense === :max ? -1 : 1
    for vals in Iterators.product(ranges...)
        l, u = copy(lb), copy(ub)
        for (k, j) in enumerate(ints)
            l[j] = u[j] = vals[k]
        end
        sol = solve(Simplexact.with_bounds(lp, l, u))
        sol.status === :optimal || continue
        (best === nothing || s * sol.objective < s * best) && (best = sol.objective)
    end
    best === nothing ? (:infeasible, nothing) : (:optimal, best)
end

@testset "knapsack" begin
    values = [60, 100, 120, 80, 30, 70]
    weights = [10, 20, 30, 25, 5, 15]
    lp, xs = knapsack(values, weights, 55)
    relax = solve(lp)
    @test relax.status === :optimal && !all(Simplexact.isint, relax.x)
    s = solve_mip(lp)
    @test s.status === :optimal
    @test s.objective == knapsack_brute(values, weights, 55)
    @test sum(weights[i] * s[xs[i]] for i in eachindex(xs)) <= 55
    @test all(Simplexact.isint, s.x)
    @test verify_mip(lp, s)
    @test s.objective <= relax.objective
    for strategy in (:best_first, :depth_first), branching in (:most_fractional, :first)
        s2 = solve_mip(lp; strategy, branching)
        @test s2.objective == s.objective && verify_mip(lp, s2)
    end
end

@testset "objective lattice rounding prunes but stays exact" begin
    @test Simplexact.objective_step(knapsack([6, 10, 4], [1, 2, 3], 4)[1]) == 2
    lp = LP(:min)
    x = add_variable!(lp, "x"; integer = true)
    y = add_variable!(lp, "y")
    set_objective!(lp, [x => 1, y => 1])
    @test Simplexact.objective_step(lp) === nothing      # y is continuous
    lp2 = LP(:min)
    a = add_variable!(lp2, "a"; integer = true)
    b = add_variable!(lp2, "b"; integer = true)
    set_objective!(lp2, [a => 1 // 2, b => 3 // 4])
    @test Simplexact.objective_step(lp2) == 1 // 4
    @test Simplexact.lattice_ceil(lp2, Q(1, 4), Q(1, 3)) == 1 // 2
end

@testset "integer infeasibility and other statuses" begin
    # 2x = 1 has a rational but no integer solution
    lp = LP()
    x = add_variable!(lp, "x"; integer = true, lb = -10, ub = 10)
    add_constraint!(lp, [x => 2], :eq, 1)
    @test solve(lp).status === :optimal
    s = solve_mip(lp)
    @test s.status === :infeasible
    @test verify_mip(lp, s)
    # parity: x + y = 2z + 1 and x = y has no integer solution
    lp2 = LP()
    v = [add_variable!(lp2, n; integer = true, ub = 20) for n in ("x", "y", "z")]
    add_constraint!(lp2, [v[1] => 1, v[2] => 1, v[3] => -2], :eq, 1)
    add_constraint!(lp2, [v[1] => 1, v[2] => -1], :eq, 0)
    s2 = solve_mip(lp2)
    @test s2.status === :infeasible && verify_mip(lp2, s2)
    # bounds that round to nothing
    lp3 = LP()
    add_variable!(lp3, "x"; lb = 1 // 3, ub = 2 // 3, integer = true)
    @test solve_mip(lp3).status === :infeasible
    # unbounded relaxation is reported honestly
    lp4 = LP(:max)
    y = add_variable!(lp4, "y"; integer = true)
    set_objective!(lp4, [y => 1])
    @test solve_mip(lp4).status === :relaxation_unbounded
    # node limit keeps the best incumbent and a valid bound
    values = [23, 31, 29, 44, 53, 38, 63, 85, 89, 82]
    weights = [92, 57, 49, 68, 60, 43, 67, 84, 87, 72]
    kp, _ = knapsack(values, weights, 165)
    lim = solve_mip(kp; node_limit = 3)
    @test lim.status === :node_limit
    full = solve_mip(kp)
    @test full.objective == knapsack_brute(values, weights, 165)
    lim.objective === nothing || @test lim.objective <= full.objective
    @test lim.bound >= full.objective
    @test !verify_mip(kp, lim)
    @test_throws ArgumentError solve_mip(kp; strategy = :random)
end

@testset "tampered proofs are rejected" begin
    lp, _ = knapsack([60, 100, 120], [10, 20, 30], 50)
    s = solve_mip(lp)
    @test s.objective == 220 && verify_mip(lp, s)
    worse = MIPSolution(:optimal, Q(160), Q[1, 0, 1], s.bound, s.nodes, s.lp_iterations, s.proof)
    @test !verify_mip(lp, worse)               # 160 is feasible but the tree proves better exists
    fake = MIPSolution(:optimal, Q(280), Q[1, 1, 1], s.bound, s.nodes, s.lp_iterations, s.proof)
    @test !verify_mip(lp, fake)                # infeasible incumbent
    lone = MIPSolution(:optimal, s.objective, s.x, s.bound, s.nodes, s.lp_iterations, ProofNode(:cutoff, solve(lp)))
    @test !verify_mip(lp, lone)                # root LP bound alone does not prove 220
    noproof = MIPSolution(:optimal, s.objective, s.x, s.bound, 0, 0, nothing)
    @test !verify_mip(lp, noproof)
end

@testset "random MIPs agree with exhaustive enumeration" begin
    rng = MersenneTwister(5)
    seen = Dict(:optimal => 0, :infeasible => 0)
    branched = 0
    for t in 1:300
        lp = random_lp(rng; boxed = true, maxn = 4, maxm = 3)
        nint = 0
        for j in 1:nvars(lp)
            if rand(rng) < 0.7
                lp.integer[j] = true
                nint += 1
            end
        end
        # use fractional row data so the LP optimum is rarely integral
        for (i, r) in enumerate(lp.rows)
            lp.rows[i] = Simplexact.Row(r.name, r.coeffs, r.sense, r.rhs + Q(rand(rng, 1:4), 5))
        end
        s = solve_mip(lp)
        status, best = mip_brute(lp)
        @test s.status === status
        status === :optimal && @test s.objective == best
        reason = check_mip_certificate(lp, s)
        @test reason === nothing
        seen[s.status] += 1
        s.nodes > 1 && (branched += 1)
    end
    @test seen[:optimal] > 80 && seen[:infeasible] > 20
    @test branched >= 20
end
