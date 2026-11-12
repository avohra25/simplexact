const TESTPROB = """
* the classic example from the MPS format documentation
NAME          TESTPROB
ROWS
 N  COST
 L  LIM1
 G  LIM2
 E  MYEQN
COLUMNS
    XONE      COST         1   LIM1         1
    XONE      LIM2         1
    YTWO      COST         2   LIM1         1
    YTWO      MYEQN       -1
    ZTHREE    COST         3   LIM2         1
    ZTHREE    MYEQN        1
RHS
    RHS1      LIM1         4   LIM2         1
    RHS1      MYEQN        7
BOUNDS
 UP BND1      XONE         4
 LO BND1      YTWO        -1
 UP BND1      YTWO         1
ENDATA
"""

# one line of fixed MPS: field 1 in columns 2-3, then 5-12, 15-22, 25-36, 40-47, 50-61
fx(f1, f2 = "", f3 = "", f4 = "", f5 = "", f6 = "") =
    rstrip(" " * rpad(f1, 2) * " " * rpad(f2, 8) * "  " * rpad(f3, 8) * "  " * rpad(f4, 12) * "   " * rpad(f5, 8) * "  " * f6)

@testset "MPS parsing" begin
    lp = parse_mps(TESTPROB)
    @test lp.sense === :min && lp.names == ["XONE", "YTWO", "ZTHREE"]
    @test lp.obj == [1, 2, 3] && iszero(lp.objconst)
    @test [r.name for r in lp.rows] == ["LIM1", "LIM2", "MYEQN"]
    @test [r.sense for r in lp.rows] == [:le, :ge, :eq] && [r.rhs for r in lp.rows] == [4, 1, 7]
    @test lp.rows[3].coeffs == Dict(2 => Q(-1), 3 => Q(1))
    @test lp.lb == [0, -1, 0] && lp.ub == [4, 1, nothing]
    s = solve(lp)
    @test s.status === :optimal && s.objective == 16 && s.x == [0, -1, 6] && verify(lp, s)

    # OBJSENSE, objective constant, free rows, exact decimals, both marker spellings
    lp = parse_mps("""
    NAME demo
    OBJSENSE
        MAX
    ROWS
     N  profit
     N  ignored
     L  cap
    COLUMNS
        MARKER  'MARKER'  'INTORG'
        x  profit  0.1  cap  1
        x  ignored  99
        MARKER  'MARKER'  'INTEND'
        y  profit  1/3  cap  1
    RHS
        rhs  profit  -2.5  cap  7
        rhs  ignored  4
    ENDATA
    """)
    @test lp.sense === :max && lp.obj == [1 // 10, 1 // 3] && lp.objconst == 5 // 2
    @test nrows(lp) == 1 && lp.integer == [true, false] && lp.rows[1].rhs == 7
    @test parse_mps("NAME\nOBJSENSE MAXIMIZE\nROWS\n N o\nCOLUMNS\n x o 1\nENDATA\n").sense === :max

    # every bound type, including the negative-UP convention and infinities
    lp = parse_mps("""
    NAME
    ROWS
     N  o
    COLUMNS
        a o 1
        b o 1
        c o 1
        d o 1
        e o 1
        f o 1
        g o 1
        h o 1
        i o 1
        k o 1
    BOUNDS
     UP BND a -3
     LO BND b -1
     UP BND b -1/2
     FX BND c 2.25
     FR BND d
     MI BND e
     UP BND e 4
     PL BND f
     BV BND g
     LI BND h -2
     UI BND h 9
     LO BND i -1e30
     UP BND i Infinity
     UP k 6
    ENDATA
    """)
    @test lp.lb == [nothing, -1, 9 // 4, nothing, nothing, 0, 0, -2, nothing, 0]
    @test lp.ub == [-3, -1 // 2, 9 // 4, nothing, 4, nothing, 1, 9, nothing, 6]
    @test lp.integer == [false, false, false, false, false, false, true, true, false, false]

    # ranges on L, G and E rows (positive and negative)
    lp = parse_mps("""
    NAME
    ROWS
     N o
     L l
     G g
     E ep
     E en
    COLUMNS
        x o 1 l 1
        x g 1 ep 1
        x en 1
    RHS
        rhs l 5 g 1
        rhs ep 2 en 3
    RANGES
        rng l 2 g -3
        rng ep 4 en -1
    ENDATA
    """)
    sig = [(r.name, r.sense, r.rhs) for r in lp.rows]
    @test sig == [("l", :le, 5), ("l_lo", :ge, 3), ("g", :ge, 1), ("g_hi", :le, 4),
                  ("ep", :ge, 2), ("ep_hi", :le, 6), ("en", :le, 3), ("en_lo", :ge, 2)]
end

@testset "fixed MPS allows spaces in names" begin
    t = join([
        "NAME          FIXED",
        "ROWS",
        fx("N", "COST"), fx("L", "LIM 1"), fx("E", "MY EQN"),
        "COLUMNS",
        fx("", "X ONE", "COST", "1", "LIM 1", "1"),
        fx("", "MARKER", "'MARKER'", "", "'INTORG'"),
        fx("", "Y TWO", "COST", "-2", "MY EQN", "1"),
        fx("", "MARKER", "'MARKER'", "", "'INTEND'"),
        "RHS",
        fx("", "", "LIM 1", "4", "MY EQN", "3"),
        "BOUNDS",
        fx("UP", "BND", "Y TWO", "5"),
        "ENDATA"], "\n")
    lp = parse_mps(t; fixed = true)
    @test lp.names == ["X ONE", "Y TWO"] && lp.integer == [false, true]
    @test [r.name for r in lp.rows] == ["LIM 1", "MY EQN"] && lp.ub == [nothing, 5]
    s = solve_mip(lp)
    @test s.status === :optimal && s.objective == -6 && verify_mip(lp, s)
    # the same text read as free MPS splits the names apart and fails
    @test_throws MPSParseError parse_mps(t)
    @test_throws ArgumentError write_mps(lp)
end

@testset "MPS errors carry line numbers" begin
    bad(text) = try
        parse_mps(text); nothing
    catch e
        e isa MPSParseError || rethrow()
        e
    end
    e = bad("NAME\nROWS\n N o\n X r\n")
    @test e.line == 4 && occursin("row type", e.msg)
    e = bad("NAME\nROWS\n N o\nCOLUMNS\n x nope 1\n")
    @test e.line == 5 && occursin("unknown row", e.msg)
    @test occursin("unknown column", bad("NAME\nROWS\n N o\nCOLUMNS\n x o 1\nBOUNDS\n UP B y 1\n").msg)
    @test occursin("SC", bad("NAME\nROWS\n N o\nCOLUMNS\n x o 1\nBOUNDS\n SC B x 1\n").msg)
    @test occursin("invalid number", bad("NAME\nROWS\n N o\nCOLUMNS\n x o abc\n").msg)
    @test occursin("only allowed in BOUNDS", bad("NAME\nROWS\n N o\nCOLUMNS\n x o 1e31\n").msg)
    @test occursin("duplicate", bad("NAME\nROWS\n N o\n L r\n G r\n").msg)
    @test occursin("objective", bad("NAME\nROWS\n L r\nCOLUMNS\n x r 1\n").msg)
    @test occursin("missing ROWS", bad("NAME x\n").msg)
    @test occursin("before ROWS", bad("NAME\nCOLUMNS\n x o 1\n").msg)
    @test occursin("unknown section", bad("NAME\nROWZ\n").msg)
    @test occursin("after ENDATA", bad("NAME\nROWS\n N o\nCOLUMNS\n x o 1\nENDATA\n x o 1\n").msg)
    @test occursin("MAX or MIN", bad("NAME\nOBJSENSE\n UP\n").msg)
    @test occursin("marker", bad("NAME\nROWS\n N o\nCOLUMNS\n M 'MARKER' 'FOO'\n").msg)
    @test bad(TESTPROB) === nothing
end

@testset "write_mps round-trips exactly" begin
    rng = MersenneTwister(2024)
    for t in 1:300
        lp = random_lp(rng; boxed = false, maxn = 4, maxm = 4)
        for j in 1:nvars(lp)
            rand(rng) < 0.3 && (lp.integer[j] = true)
            rand(rng) < 0.1 && (lp.ub[j] = Q(-rand(rng, 1:3)); lp.lb[j] = rand(rng) < 0.5 ? nothing : Q(0))
            rand(rng) < 0.1 && (lp.lb[j] = Q(0); lp.ub[j] = Q(1))
            rand(rng) < 0.1 && (lp.obj[j] = Q(rand(rng, -9:9), rand(rng, 1:8)))
        end
        rand(rng) < 0.3 && (lp.objconst = Q(rand(rng, -20:20), rand(rng, [1, 3, 4])))
        isempty(lp.rows) || (lp.rows[1] = Simplexact.Row("obj", lp.rows[1].coeffs, lp.rows[1].sense, lp.rows[1].rhs - Q(1, 7)))
        text = write_mps(lp)
        lp2 = parse_mps(text)
        @test lp2.names == lp.names && lp2.sense === lp.sense
        @test lp2.obj == lp.obj && lp2.objconst == lp.objconst
        @test lp2.lb == lp.lb && lp2.ub == lp.ub && lp2.integer == lp.integer
        @test [(r.name, r.coeffs, r.sense, r.rhs) for r in lp2.rows] == [(r.name, r.coeffs, r.sense, r.rhs) for r in lp.rows]
        @test write_mps(lp2) == text
        # and MPS and LP format describe the same problem
        @test write_lp(parse_lp(write_lp(lp2))) == write_lp(parse_lp(write_lp(lp)))
    end
end

@testset "mps_value" begin
    @test Simplexact.mps_value(Q(3)) == "3"
    @test Simplexact.mps_value(Q(-1, 8)) == "-0.125"
    @test Simplexact.mps_value(Q(7, 20)) == "0.35"
    @test Simplexact.mps_value(Q(1, 3)) == "1/3"
    @test_throws ArgumentError Simplexact.mps_value(Q(1, 3); strict = true)
    for q in (Q(1, 1024), Q(-123456789, 100000), Q(5, 2), Q(-1, 5))
        @test parse_q(Simplexact.mps_value(q; strict = true)) == q
    end
end

@testset "read_model picks the format" begin
    dir = mktempdir()
    p1 = joinpath(dir, "a.mps"); write(p1, TESTPROB)
    p2 = joinpath(dir, "a.lp"); write(p2, write_lp(parse_mps(TESTPROB)))
    p3 = joinpath(dir, "noext"); write(p3, TESTPROB)
    p4 = joinpath(dir, "noext2"); write(p4, read(p2, String))
    for p in (p1, p2, p3, p4)
        @test solve(read_model(p)).objective == 16
    end
    @test read_model(p3; format = :mps).names == ["XONE", "YTWO", "ZTHREE"]
    @test_throws ArgumentError read_model(p1; format = :xls)
end
