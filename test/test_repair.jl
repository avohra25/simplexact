@testset "feasibility repair on the staffing plan" begin
    lp = read_lp(joinpath(EXAMPLES, "staffing.lp"))
    r = feasibility_repair(lp)
    @test r.status === :optimal && r.violation == 3 // 4
    @test r.row_violation == [3 // 4, 0, 0, 0, 0, 0]       # only coverage gives way
    @test all(iszero, r.lower_violation) && all(iszero, r.upper_violation)
    @test verify_repair(lp, r)
    # make coverage expensive and mentoring gives way; make both expensive and the budget does
    r = feasibility_repair(lp; weights = Dict("coverage" => 1000))
    @test r.row_violation == [0, 3, 0, 0, 0, 0] && r.violation == 3 && verify_repair(lp, r)
    r = feasibility_repair(lp; weights = Dict("coverage" => 1000, "mentoring" => 1000))
    @test r.row_violation[1] == 0 && r.row_violation[3] == 60 && r.violation == 60
    @test verify_repair(lp, r)
    @test feasibility_repair(lp; weights = fill(2, 6)).violation == 3 // 2
    @test_throws ArgumentError feasibility_repair(lp; weights = Dict("nope" => 1))
    @test_throws ArgumentError feasibility_repair(lp; weights = [1, 2])
    @test_throws ArgumentError feasibility_repair(lp; weights = fill(0, 6))
    @test_throws ArgumentError feasibility_repair(lp; relax_bounds = true, bound_weight = 0)
end

@testset "repair with relaxed bounds" begin
    lp = LP(:max)
    x = add_variable!(lp, "x"; ub = 2)
    add_constraint!(lp, [x => 1], :ge, 5; name = "need")
    r = feasibility_repair(lp)
    @test r.violation == 3 && r.x == [2] && verify_repair(lp, r)
    r = feasibility_repair(lp; relax_bounds = true, bound_weight = "1/2")
    @test r.violation == 3 // 2 && r.x == [5] && r.upper_violation == [3] && r.row_violation == [0]
    @test verify_repair(lp, r)
    # contradictory bounds: hopeless unless bounds may move
    add_variable!(lp, "y"; lb = 4, ub = 1)
    @test feasibility_repair(lp).status === :infeasible
    r = feasibility_repair(lp; relax_bounds = true)
    @test r.status === :optimal && r.violation == 3 + 3 && verify_repair(lp, r)
end

@testset "optimize among least-violation points" begin
    lp = LP(:max)
    x = add_variable!(lp, "x")
    y = add_variable!(lp, "y")
    set_objective!(lp, [x => 1, y => 2])
    add_constraint!(lp, [x => 1, y => 1], :le, 4; name = "cap")
    add_constraint!(lp, [x => 1, y => 1], :ge, 6; name = "floor")
    add_constraint!(lp, [y => 1], :le, 3; name = "ycap")
    r = feasibility_repair(lp; optimize = true)
    @test r.violation == 2 && r.objective == 9 && r.x == [3, 3]     # anywhere on 4 <= x+y <= 6 costs 2
    @test verify_repair(lp, r)
    # unbounded among the least-violation points
    lp2 = LP(:max)
    a = add_variable!(lp2, "a")
    b = add_variable!(lp2, "b")
    set_objective!(lp2, [a => 1])
    add_constraint!(lp2, [b => 1], :le, -1)
    r = feasibility_repair(lp2; optimize = true)
    @test r.status === :unbounded && r.violation == 1
    # feasible problems repair to their optimum
    for f in ("diet.lp", "transport.lp", "knapsack.lp")
        p = read_lp(joinpath(EXAMPLES, f))
        r = feasibility_repair(p; optimize = true)
        @test r.violation == 0 && verify_repair(p, r)
        @test r.objective == (any(p.integer) ? solve_mip(p) : solve(p)).objective
    end
end

@testset "integer variables stay integer" begin
    lp = LP()
    x = add_variable!(lp, "x"; integer = true, ub = 10)
    add_constraint!(lp, [x => 2], :eq, 3; name = "odd")
    r = feasibility_repair(lp)
    @test r.violation == 1 && isone(denominator(r.x[1])) && verify_repair(lp, r)
    @test feasibility_repair(LP()).violation == 0
end

@testset "random repairs are least-violation and verify" begin
    rng = MersenneTwister(77)
    for t in 1:400
        lp = random_lp(rng; boxed = true, maxn = 3, maxm = 4)
        ints = rand(rng) < 0.25
        ints && (lp.integer .= rand(rng, Bool, nvars(lp)))
        w = Q.(rand(rng, 1:4, nrows(lp)))
        relax = rand(rng) < 0.3
        r = feasibility_repair(lp; weights = w, relax_bounds = relax, optimize = rand(rng, Bool))
        r.status === :optimal || continue
        @test verify_repair(lp, r)
        feas = (ints ? solve_mip(lp) : solve(lp)).status === :optimal
        @test iszero(r.violation) == feas
        if !ints && !feas && !relax
            # every IIS must have a violated member at the repaired point
            iis = find_iis(lp)
            @test any(!iszero(r.row_violation[i]) for i in iis.rows) ||
                  any(r.x[j] < lp.lb[j] for j in iis.lower) || any(r.x[j] > lp.ub[j] for j in iis.upper)
        end
        if !ints && !relax && nrows(lp) > 0 && nvars(lp) + 2 * nrows(lp) <= 7
            # brute force: no vertex of the elastic polytope beats the repair
            m = Simplexact.elastic_model(lp, w, nothing)
            for j in nvars(lp)+1:nvars(m)
                m.ub[j] = Q(sum(abs(a) for r in lp.rows for a in values(r.coeffs); init = 0) * 10 + 50)
            end
            st, best = brute_force(m)
            @test st === :optimal && best == r.violation
        end
    end
end

@testset "tampered repairs are rejected" begin
    lp = read_lp(joinpath(EXAMPLES, "staffing.lp"))
    r = feasibility_repair(lp)
    T(; kw...) = begin
        f = Dict{Symbol,Any}(n => getfield(r, n) for n in fieldnames(Repair))
        merge!(f, Dict(kw))
        Repair((f[n] for n in fieldnames(Repair))...)
    end
    @test occursin("not the proven least", check_repair(lp, T(violation = Q(1, 2))))
    @test occursin("misreported", check_repair(lp, T(row_violation = zeros(Q, 6))))
    @test occursin("does not match", check_repair(lp, T(x = [Q(0), Q(0), Q(0)])))
    @test occursin("elastic model does not match", check_repair(lp, T(weights = fill(Q(2), 6))))
    @test occursin("objective", check_repair(lp, T(objective = Q(1))))
    m = deepcopy(r.model); m.rows[1] = Simplexact.Row(m.rows[1].name, m.rows[1].coeffs, :ge, Q(11))
    @test occursin("elastic model does not match", check_repair(lp, T(model = m)))
end
