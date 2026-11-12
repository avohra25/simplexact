const EXAMPLES = joinpath(@__DIR__, "..", "examples")

@testset "LP format parsing" begin
    lp = parse_lp("""
    \\ a comment line
    Maximize
      profit: 3 x + 2y + 1/2   \\ trailing comment
    Subject To
      cap: x + y <= 4
      x + 3 * y
         =< 6
      lim: x - 3 >= -6 + 0 y
      eq: 2 x = y + 1
    Bounds
      x <= 3
      -1 <= y <= 5/2
      z free
      -inf <= w <= 7
      v >= -2.5
      u = 4
    General
      x
    Binary
      b
    End
    """)
    @test lp.sense === :max
    @test lp.names == ["x", "y", "z", "w", "v", "u", "b"]
    @test lp.obj[1:2] == [3, 2] && lp.objconst == 1 // 2
    @test nrows(lp) == 4
    @test lp.rows[1].name == "cap" && lp.rows[2].name == "c2"
    @test lp.rows[2].coeffs == Dict(1 => Q(1), 2 => Q(3)) && lp.rows[2].sense === :le && lp.rows[2].rhs == 6
    @test lp.rows[3].sense === :ge && lp.rows[3].rhs == -3 && lp.rows[3].coeffs == Dict(1 => Q(1))
    @test lp.rows[4].coeffs == Dict(1 => Q(2), 2 => Q(-1)) && lp.rows[4].rhs == 1
    @test lp.lb == [0, -1, nothing, nothing, -5 // 2, 4, 0]
    @test lp.ub == [3, 5 // 2, nothing, 7, nothing, 4, 1]
    @test lp.integer == [true, false, false, false, false, false, true]
end

@testset "LP format errors carry line numbers" begin
    bad(text) = try
        parse_lp(text); nothing
    catch e
        e isa LPParseError || rethrow()
        e
    end
    e = bad("minimize\n x\nsubject to\n x + y\n")
    @test e.line == 4 && occursin("comparison", e.msg)
    e = bad("minimize\n x\nsubject to\n -1 <= x + y <= 2\n")
    @test e.line == 4 && occursin("ranged", e.msg)
    @test occursin("section keyword", bad("x + y\n").msg)
    @test occursin("missing", bad("subject to\n x <= 1\n").msg)
    @test bad("minimize\n x ? y\n").msg == "unexpected character '?'"
    @test occursin("'+' or '-'", bad("minimize\n x y\n").msg)
    @test occursin("bound", bad("minimize\n x\nbounds\n x <= y\n").msg)
    @test occursin("after 'end'", bad("minimize\n x\nend\n x\n").msg)
    @test occursin("duplicate", bad("min\n x\nst\n r: x <= 1\n r: x >= 0\n").msg)
    @test occursin("infinity", bad("min\n x\nst\n x + inf <= 2\n").msg)
    @test occursin("zero denominator", bad("min\n x\nst\n r: 2 x <= 1/0\n").msg)
    @test bad("max\n x\nmin\n x\n").line == 3
    @test bad("max\n x\nst\n x <= 1\nend\n") === nothing
end

@testset "write_lp round-trips exactly" begin
    rng = MersenneTwister(99)
    for t in 1:200
        lp = random_lp(rng; boxed = false, maxn = 4, maxm = 4)
        for j in 1:nvars(lp)
            rand(rng) < 0.3 && (lp.integer[j] = true)
        end
        isempty(lp.rows) || (lp.rows[1] = Simplexact.Row("first", lp.rows[1].coeffs, lp.rows[1].sense, lp.rows[1].rhs - Q(1, 7)))
        lp2 = parse_lp(write_lp(lp))
        # the parser numbers variables by first appearance, so compare by name
        perm = [findfirst(==(n), lp2.names) for n in lp.names]
        @test sort(lp2.names) == sort(lp.names) && lp2.sense === lp.sense
        @test lp2.obj[perm] == lp.obj && lp2.objconst == lp.objconst
        @test lp2.lb[perm] == lp.lb && lp2.ub[perm] == lp.ub && lp2.integer[perm] == lp.integer
        rowsig(m, r) = (r.name, Dict(m.names[j] => a for (j, a) in r.coeffs), r.sense, r.rhs)
        @test [rowsig(lp2, r) for r in lp2.rows] == [rowsig(lp, r) for r in lp.rows]
        @test write_lp(parse_lp(write_lp(lp2))) == write_lp(lp2)
        a, b = solve(lp), solve(lp2)
        @test a.status === b.status && a.objective == b.objective
    end
end

@testset "examples" begin
    diet = read_lp(joinpath(EXAMPLES, "diet.lp"))
    s = solve(diet)
    @test s.status === :optimal && s.objective == 25603 // 97 && verify(diet, s)
    beale = read_lp(joinpath(EXAMPLES, "beale.lp"))
    @test solve(beale).objective == -5 // 4
    tr = read_lp(joinpath(EXAMPLES, "transport.lp"))
    st = solve(tr)
    @test st.objective == 3415 // 2 && verify(tr, st)
    inf = read_lp(joinpath(EXAMPLES, "infeasible.lp"))
    si = solve(inf)
    @test si.status === :infeasible && verify(inf, si)
    unb = read_lp(joinpath(EXAMPLES, "unbounded.lp"))
    su = solve(unb)
    @test su.status === :unbounded && verify(unb, su)
    kp = read_lp(joinpath(EXAMPLES, "knapsack.lp"))
    sk = solve_mip(kp)
    @test sk.objective == knapsack_brute([23, 31, 29, 44, 53, 38, 63, 85, 89, 82],
        [92, 57, 49, 68, 60, 43, 67, 84, 87, 72], 165)
    @test verify_mip(kp, sk)
end

@testset "CLI" begin
    run_cli(args...) = begin
        out, err = IOBuffer(), IOBuffer()
        code = Simplexact.main(collect(String, args); io = out, err = err)
        code, String(take!(out)), String(take!(err))
    end
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "diet.lp"), "--decimal", "--duals")
    @test code == 0
    @test occursin("status: optimal", out) && occursin("objective: 25603/97", out)
    @test occursin("≈ 263.948454", out) && occursin("shadow prices:", out)
    @test occursin("certificate: verified", out)
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "knapsack.lp"), "--strategy", "depth_first")
    @test code == 0 && occursin("proof: verified", out)
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "knapsack.lp"), "--relax")
    @test code == 0 && occursin("certificate: verified", out) && !occursin("nodes:", out)
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "knapsack.lp"), "--node-limit", "2")
    @test code == 0 && occursin("status: node_limit", out) && occursin("bound:", out)
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "infeasible.lp"))
    @test code == 0 && occursin("Farkas", out)
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "unbounded.lp"))
    @test code == 0 && occursin("improving ray", out)
    code, out, _ = run_cli("check", joinpath(EXAMPLES, "knapsack.lp"))
    @test code == 0 && occursin("10 variables (10 integer), 1 rows", out)
    code, out, _ = run_cli("format", joinpath(EXAMPLES, "transport.lp"))
    @test code == 0 && parse_lp(out).names == read_lp(joinpath(EXAMPLES, "transport.lp")).names
    path = tempname() * ".lp"
    write(path, "minimize\n x\nsubject to\n x +\n")
    code, _, err = run_cli("solve", path)
    @test code == 1 && occursin("line 4", err)
    @test run_cli()[1] == 1
    @test run_cli("solve", "/nonexistent.lp")[1] == 1
    @test run_cli("solve", joinpath(EXAMPLES, "diet.lp"), "--rule", "nope")[1] == 1
    @test run_cli("solve", joinpath(EXAMPLES, "diet.lp"), "--bogus")[1] == 1
    @test run_cli("help")[1] == 0
    # MPS input and output
    code, out, _ = run_cli("solve", joinpath(EXAMPLES, "testprob.mps"))
    @test code == 0 && occursin("objective: 37/2", out) && occursin("certificate: verified", out)
    code, out, _ = run_cli("check", joinpath(EXAMPLES, "testprob.mps"))
    @test code == 0 && occursin("3 variables (0 integer), 4 rows", out)
    code, out, _ = run_cli("format", joinpath(EXAMPLES, "knapsack.lp"), "--to", "mps")
    @test code == 0 && startswith(out, "NAME          KNAPSACK") && occursin("'INTORG'", out)
    @test solve_mip(parse_mps(out)).objective == solve_mip(read_lp(joinpath(EXAMPLES, "knapsack.lp"))).objective
    code, out, _ = run_cli("format", joinpath(EXAMPLES, "testprob.mps"))
    @test code == 0 && occursin("LIM1_lo: XONE + YTWO >= 3/2", out)
    mps = tempname() * ".txt"
    write(mps, read(joinpath(EXAMPLES, "testprob.mps"), String))
    @test run_cli("solve", mps, "--format", "mps")[1] == 0
    code, _, err = run_cli("solve", mps, "--format", "lp")
    @test code == 1 && occursin("line", err)
    @test run_cli("solve", mps, "--format", "xml")[1] == 1
    thirds = tempname() * ".lp"
    write(thirds, "min\n 1/3 x\nst\n x >= 1/3\nend\n")
    @test run_cli("format", thirds, "--to", "mps")[1] == 0
    code, _, err = run_cli("format", thirds, "--to", "mps", "--strict")
    @test code == 1 && occursin("exact decimal", err)
    @test run_cli("format", thirds, "--to", "json")[1] == 1
end
