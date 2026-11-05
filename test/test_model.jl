@testset "numbers" begin
    @test parse_q("3") == 3
    @test parse_q("-7/2") == -7 // 2
    @test parse_q("0.125") == 1 // 8
    @test parse_q("1.5e-3") == 3 // 2000
    @test parse_q("+2e2") == 200
    @test parse_q(".5") == 1 // 2
    @test_throws ArgumentError parse_q("abc")
    @test_throws ArgumentError parse_q("1/0")
    @test_throws ArgumentError parse_q(".")
    @test fmt(Q(6, 4)) == "3/2"
    @test fmt(Q(-4, 2)) == "-2"
    @test decimal(Q(1, 3)) == "0.333333"
    @test decimal(Q(-2, 3); digits = 2) == "-0.67"
    @test decimal(Q(5)) == "5.000000"
    @test decimal(Q(1, 1000); digits = 2) == "0.00"
end

@testset "model building" begin
    lp = LP(:max)
    x = add_variable!(lp, "x")
    y = add_variable!(lp, "y"; ub = 5, lb = -1)
    @test nvars(lp) == 2
    @test lp["y"] == y
    @test_throws ArgumentError add_variable!(lp, "x")
    @test_throws ArgumentError LP(:best)
    set_objective!(lp, [x => 3, y => "1/2", x => 1]; constant = 2)
    @test lp.obj == [4, 1 // 2]
    @test lp.objconst == 2
    i = add_constraint!(lp, [x => 1, y => 1], "<=", 4; name = "cap")
    @test lp.rows[i].name == "cap" && lp.rows[i].sense === :le
    @test_throws ArgumentError add_constraint!(lp, [x => 1], :lt, 1)
    @test_throws ArgumentError add_constraint!(lp, [Var(9) => 1], :le, 1)
    @test is_feasible(lp, Q[1, 2])
    @test !is_feasible(lp, Q[3, 2])
    @test !is_feasible(lp, Q[1, 6])
    @test !is_feasible(lp, Q[-1, 0])
    @test objective_value(lp, Q[1, 2]) == 4 + 1 + 2
    @test tonum(0.1) == 1 // 10
end
