function lp_of(sense, vars, obj, rows)
    lp = LP(sense)
    vs = Var[]
    for (name, lb, ub) in vars
        push!(vs, add_variable!(lp, name; lb, ub))
    end
    set_objective!(lp, [vs[j] => c for (j, c) in enumerate(obj)])
    for (coeffs, s, rhs) in rows
        add_constraint!(lp, [vs[j] => a for (j, a) in enumerate(coeffs)], s, rhs)
    end
    lp
end

@testset "textbook optimum" begin
    # max 3x + 2y  s.t. x + y <= 4, x + 3y <= 6, x <= 3
    lp = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [3, 2],
        [([1, 1], :le, 4), ([1, 3], :le, 6), ([1, 0], :le, 3)])
    s = solve(lp)
    @test s.status === :optimal
    @test s.objective == 11
    @test s.x == [3, 1]
    @test is_feasible(lp, s.x)
end

@testset "fractional optimum is exact" begin
    # max x + y s.t. 2x + y <= 3, x + 3y <= 4  ->  x = 1, y = 1 ... use thirds
    lp = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [1, 1],
        [([3, 1], :le, 2), ([1, 3], :le, 1)])
    s = solve(lp)
    @test s.status === :optimal
    @test s.x == [5 // 8, 1 // 8]
    @test s.objective == 3 // 4
end

@testset "equalities, >= rows and negative right-hand sides" begin
    # min 2x + 3y s.t. x + y >= 10, x - y = 2, y >= 1
    lp = lp_of(:min, [("x", 0, nothing), ("y", 1, nothing)], [2, 3],
        [([1, 1], :ge, 10), ([1, -1], :eq, 2)])
    s = solve(lp)
    @test s.status === :optimal
    @test s.x == [6, 4]
    @test s.objective == 24
    # same thing with a negated row
    lp2 = lp_of(:min, [("x", 0, nothing), ("y", 1, nothing)], [2, 3],
        [([-1, -1], :le, -10), ([-1, 1], :eq, -2)])
    @test solve(lp2).objective == 24
end

@testset "bounds: negative lower, upper only, free" begin
    lp = lp_of(:min, [("x", -5, 5), ("y", nothing, 3), ("z", nothing, nothing)], [1, -1, 2],
        [([1, 1, 1], :ge, -2), ([0, 0, 1], :le, 10), ([0, 0, 1], :ge, -10)])
    s = solve(lp)
    @test s.status === :optimal
    @test is_feasible(lp, s.x)
    # z = -10 is worth -20, which forces x + y >= 8, so x = 5 and y = 3
    @test s.x == [5, 3, -10]
    @test s.objective == -18
end

@testset "statuses" begin
    infeasible = lp_of(:min, [("x", 0, nothing)], [1], [([1], :le, 1), ([1], :ge, 2)])
    @test solve(infeasible).status === :infeasible
    clash = lp_of(:min, [("x", 3, 2)], [1], [])
    @test solve(clash).status === :infeasible
    unbounded = lp_of(:max, [("x", 0, nothing), ("y", 0, nothing)], [1, 1], [([1, -1], :le, 1)])
    s = solve(unbounded)
    @test s.status === :unbounded
    @test s.ray[1] + s.ray[2] > 0
    @test s.ray[1] - s.ray[2] <= 0
    # free variable with no constraints
    free = lp_of(:min, [("x", nothing, nothing)], [1], [])
    @test solve(free).status === :unbounded
    # no rows at all, bounded by variable bounds
    box = lp_of(:max, [("x", 0, 4), ("y", -2, 3)], [2, -1], [])
    sb = solve(box)
    @test sb.status === :optimal && sb.x == [4, -2] && sb.objective == 10
    empty = LP()
    @test solve(empty).status === :optimal
end

@testset "redundant and degenerate" begin
    # duplicated equality rows
    lp = lp_of(:min, [("x", 0, nothing), ("y", 0, nothing)], [1, 2],
        [([1, 1], :eq, 4), ([2, 2], :eq, 8), ([1, 1], :eq, 4)])
    s = solve(lp)
    @test s.status === :optimal && s.objective == 4 && s.x == [4, 0]
    # inconsistent duplicates
    bad = lp_of(:min, [("x", 0, nothing), ("y", 0, nothing)], [1, 2],
        [([1, 1], :eq, 4), ([2, 2], :eq, 9)])
    @test solve(bad).status === :infeasible
    # Beale's cycling example: Dantzig's rule alone can loop, Bland and :auto must not
    beale = lp_of(:min, [(n, 0, nothing) for n in ("a", "b", "c", "d")], [-3 // 4, 20, -1 // 2, 6],
        [([1 // 4, -8, -1, 9], :le, 0), ([1 // 2, -12, -1 // 2, 3], :le, 0), ([0, 0, 1, 0], :le, 1)])
    for rule in (:auto, :bland)
        s = solve(beale; rule)
        @test s.status === :optimal
        @test s.objective == -5 // 4
        @test is_feasible(beale, s.x)
    end
    @test_throws ArgumentError solve(beale; rule = :nope)
end
