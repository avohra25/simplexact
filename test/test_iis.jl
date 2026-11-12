@testset "IIS of the textbook infeasible LP" begin
    lp = read_lp(joinpath(EXAMPLES, "infeasible.lp"))
    iis = find_iis(lp)
    @test iis.rows == [1, 2, 3] && isempty(iis.lower) && isempty(iis.upper)
    @test iis.farkas == [2 // 5, 1 // 5, -1]
    @test verify_iis(lp, iis)
end

@testset "IIS picks out the conflicting part" begin
    lp = LP(:max)
    x = add_variable!(lp, "x"; ub = 10)
    y = add_variable!(lp, "y")
    z = add_variable!(lp, "z"; lb = -5, ub = 5)
    set_objective!(lp, [x => 1, y => 1, z => 1])
    add_constraint!(lp, [z => 1], :le, 3; name = "harmless")
    add_constraint!(lp, [x => 1, y => 1], :ge, 12; name = "need")
    add_constraint!(lp, [y => 1], :le, 1; name = "ycap")
    add_constraint!(lp, [z => 1, y => 1], :ge, -50; name = "loose")
    iis = find_iis(lp)
    @test [lp.rows[i].name for i in iis.rows] == ["need", "ycap"]
    @test iis.upper == [1] && isempty(iis.lower)            # x <= 10 is part of the reason
    @test verify_iis(lp, iis)

    # contradictory bounds are an IIS of their own
    w = add_variable!(lp, "w"; lb = 3, ub = 2)
    iis = find_iis(lp)
    @test isempty(iis.rows) && iis.lower == [w.index] && iis.upper == [w.index] && isempty(iis.farkas)
    @test verify_iis(lp, iis)

    # equality rows and a negative lower bound
    lp = LP()
    a = add_variable!(lp, "a"; lb = -2)
    b = add_variable!(lp, "b"; lb = nothing)
    add_constraint!(lp, [a => 1, b => 1], :eq, 0)
    add_constraint!(lp, [b => 1], :ge, 3)
    add_constraint!(lp, [a => 1, b => -1], :ge, -100)
    iis = find_iis(lp)
    @test iis.rows == [1, 2] && iis.lower == [1] && verify_iis(lp, iis)
end

@testset "find_iis refuses feasible problems" begin
    @test_throws ArgumentError find_iis(read_lp(joinpath(EXAMPLES, "diet.lp")))
    @test_throws ArgumentError find_iis(read_lp(joinpath(EXAMPLES, "unbounded.lp")))
    # integer-infeasible but LP-feasible: no IIS of the relaxation
    lp = LP()
    x = add_variable!(lp, "x"; integer = true)
    add_constraint!(lp, [x => 2], :eq, 1)
    @test solve_mip(lp).status === :infeasible
    @test_throws ArgumentError find_iis(lp)
end

@testset "random infeasible LPs give verified, irreducible IISs" begin
    rng = MersenneTwister(4242)
    found = 0
    for t in 1:1200
        lp = random_lp(rng; boxed = rand(rng, Bool), maxn = 4, maxm = 5)
        solve(lp).status === :infeasible || continue
        found += 1
        iis = find_iis(lp)
        @test verify_iis(lp, iis)
        # irreducibility, checked by re-solving rather than trusting the witnesses
        rows, lo, up = iis.rows, Set(iis.lower), Set(iis.upper)
        sub(r, l, u) = Simplexact.subsystem(lp, r, l, u)
        @test solve(sub(rows, lo, up)).status === :infeasible
        for i in rows
            @test solve(sub(filter(!=(i), rows), lo, up)).status === :optimal
        end
        for j in lo
            @test solve(sub(rows, setdiff(lo, (j,)), up)).status === :optimal
        end
        for j in up
            @test solve(sub(rows, lo, setdiff(up, (j,)))).status === :optimal
        end
        # all pivoting rules find an IIS that verifies (not necessarily the same one)
        @test verify_iis(lp, find_iis(lp; rule = :bland))
    end
    @test found > 100
end

@testset "tampered IISs are rejected" begin
    lp = LP(:min)
    x = add_variable!(lp, "x")
    y = add_variable!(lp, "y")
    add_constraint!(lp, [x => 1, y => 1], :le, 1; name = "small")
    add_constraint!(lp, [x => 1], :ge, 2; name = "bigx")
    add_constraint!(lp, [y => 1], :le, 7; name = "spare")
    iis = find_iis(lp)
    @test iis.rows == [1, 2] && iis.lower == [2] && verify_iis(lp, iis)
    T(; rows = iis.rows, lower = iis.lower, upper = iis.upper, farkas = iis.farkas, witnesses = iis.witnesses) =
        Simplexact.IIS(rows, lower, upper, farkas, witnesses)
    @test occursin("not proven infeasible", check_iis(lp, T(rows = [1], farkas = iis.farkas[1:1])))
    @test occursin("not proven infeasible", check_iis(lp, T(lower = Int[])))
    @test occursin("not proven infeasible", check_iis(lp, T(farkas = -iis.farkas)))
    # an extra member that is not needed has no witness
    @test occursin("no witness", check_iis(lp, T(rows = [1, 2, 3], farkas = vcat(iis.farkas, 0))))
    bad = copy(iis.witnesses); bad[(:row, 1)] = [Q(0), Q(0)]
    @test occursin("violates", check_iis(lp, T(witnesses = bad)))
    @test occursin("no finite upper bound", check_iis(lp, T(upper = [1])))
    @test occursin("out of range", check_iis(lp, T(rows = [1, 9])))
    @test occursin("repeated", check_iis(lp, T(rows = [1, 1])))
end
