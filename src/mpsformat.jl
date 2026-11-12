# Reading and writing MPS files, the column-oriented format most LP and MIP
# solvers exchange problems in.
#
#   NAME          TESTPROB
#   OBJSENSE                     (optional; MAX or MIN on this line or the next)
#       MAX
#   ROWS
#    N  COST
#    L  LIM1
#    G  LIM2
#    E  MYEQN
#   COLUMNS
#       MARKER    'MARKER'   'INTORG'
#       XONE      COST  1    LIM1  1
#       MARKER    'MARKER'   'INTEND'
#       YTWO      COST  2    MYEQN -1
#   RHS
#       RHS1      LIM1  4    MYEQN 7
#   RANGES
#       RNG1      LIM1  2.5
#   BOUNDS
#    UP BND1      XONE  4
#    MI BND1      YTWO
#   ENDATA
#
# Both free MPS (fields separated by whitespace) and fixed MPS (fields in
# fixed columns, so names may contain spaces) are read. Numbers are parsed
# exactly, so "0.1" is 1/10.

"Error raised for malformed MPS text, with the 1-based line number."
struct MPSParseError <: Exception
    line::Int
    msg::String
end
Base.showerror(io::IO, e::MPSParseError) = print(io, "line ", e.line, ": ", e.msg)

const MPS_SECTIONS = ("NAME", "OBJSENSE", "OBJSENCE", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "ENDATA")

# Fixed MPS field columns (1-based, inclusive): field 1 is the row type or
# bound type, fields 2-6 are name, name, number, name, number.
const FIXED_FIELDS = ((2, 3), (5, 12), (15, 22), (25, 36), (40, 47), (50, 61))

function fixed_fields(line::AbstractString)
    s = rstrip(line)
    chars = collect(s)
    out = String[]
    for (a, b) in FIXED_FIELDS
        a > length(chars) && break
        push!(out, strip(String(chars[a:min(b, end)])))
    end
    while !isempty(out) && isempty(out[end])
        pop!(out)
    end
    out
end

"Values this large (or the words Inf/Infinity) are treated as infinite bounds, as in most solvers."
const MPS_INFINITY = Q(10)^30

function mps_number(s::AbstractString, ln::Int)
    t = lowercase(s)
    t in ("inf", "+inf", "infinity", "+infinity") && return :posinf
    t in ("-inf", "-infinity") && return :neginf
    q = try
        parse_q(s)
    catch e
        e isa ArgumentError || rethrow()
        throw(MPSParseError(ln, e.msg))
    end
    q >= MPS_INFINITY ? :posinf : q <= -MPS_INFINITY ? :neginf : q
end

function finite_number(s::AbstractString, ln::Int)
    v = mps_number(s, ln)
    v isa Symbol && throw(MPSParseError(ln, "infinite value \"$s\" is only allowed in BOUNDS"))
    v
end

"""
    parse_mps(text; fixed = false) -> LP

Parse MPS text. With `fixed = false` (the default) fields are separated by
whitespace, which reads free MPS and any fixed MPS file whose names contain no
spaces. With `fixed = true` fields are taken from the classic fixed columns.

Conventions, matching CPLEX and Gurobi:

* the first `N` row is the objective; further `N` rows are free rows and are
  dropped;
* a right-hand side on the objective row is minus the objective constant;
* variables default to `0 <= x < inf`; an `UP` bound below zero on a variable
  whose lower bound was not set makes the lower bound `-inf`;
* integer variables between `INTORG`/`INTEND` markers keep the default
  bounds `[0, inf)`;
* a `RANGES` entry turns a row into two rows: the original keeps its name and
  one side, and a second row named `<row>_lo` or `<row>_hi` holds the other.
"""
function parse_mps(text::AbstractString; fixed::Bool = false)
    lp = LP(:min)
    section = :none
    objrow = ""
    freerows = Set{String}()
    rowindex = Dict{String,Int}()
    colindex = Dict{String,Var}()
    inint = false
    rhsset = nothing
    rngset = nothing
    bndset = nothing
    lbset = Set{Int}()
    ranges = Dict{Int,Tuple{Q,Int}}()          # row => (range value, line)
    sawrows = false
    ended = false
    rows_coeffs = Pair{String,Vector{Pair{Int,Q}}}[]

    for (ln, raw) in enumerate(split(text, '\n'))
        line = replace(String(raw), '\r' => "")
        isempty(strip(line)) && continue
        startswith(line, '*') && continue
        ended && throw(MPSParseError(ln, "text after ENDATA"))
        if !isspace(line[1])
            toks = split(line)
            key = uppercase(toks[1])
            key in MPS_SECTIONS || throw(MPSParseError(ln, "unknown section \"$(toks[1])\""))
            if key == "NAME"
                section = :name
            elseif key in ("OBJSENSE", "OBJSENCE")
                section = :objsense
                if length(toks) >= 2
                    set_mps_sense!(lp, toks[2], ln)
                    section = :none
                end
            elseif key == "ROWS"
                sawrows = true
                section = :rows
            elseif key == "ENDATA"
                ended = true
            else
                sawrows || throw(MPSParseError(ln, "$key before ROWS"))
                section = Symbol(lowercase(key))
            end
            length(toks) > 2 && key != "NAME" && throw(MPSParseError(ln, "unexpected text after $key"))
            continue
        end
        f = fixed ? fixed_fields(line) : String.(split(line))
        if fixed && section in (:columns, :rhs, :ranges)
            # field 1 (columns 2-3) is blank in these sections
            (isempty(f) || isempty(f[1])) || throw(MPSParseError(ln, "unexpected text in columns 2-3"))
            f = f[2:end]
        end
        isempty(f) && continue
        if section === :objsense
            set_mps_sense!(lp, f[1], ln)
        elseif section === :rows
            length(f) == 2 || throw(MPSParseError(ln, "expected a row type and a row name"))
            t, name = uppercase(f[1]), f[2]
            (haskey(rowindex, name) || name == objrow || name in freerows) &&
                throw(MPSParseError(ln, "duplicate row name \"$name\""))
            if t == "N"
                isempty(objrow) ? (objrow = name) : push!(freerows, name)
            elseif t in ("L", "G", "E")
                add_constraint!(lp, Pair{Var,Q}[], t == "L" ? :le : t == "G" ? :ge : :eq, 0; name)
                rowindex[name] = nrows(lp)
            else
                throw(MPSParseError(ln, "unknown row type \"$(f[1])\""))
            end
        elseif section === :columns
            if length(f) >= 3 && strip(f[2], ['\'', '"']) == "MARKER"
                m = uppercase(strip(f[end], ['\'', '"']))
                m == "INTORG" ? (inint = true) : m == "INTEND" ? (inint = false) :
                    throw(MPSParseError(ln, "unknown marker \"$(f[3])\""))
                continue
            end
            (length(f) in (1, 3, 5)) || throw(MPSParseError(ln, "expected a column name and row/value pairs"))
            name = f[1]
            v = get!(() -> add_variable!(lp, name; integer = inint), colindex, name)
            inint && (lp.integer[v.index] = true)
            for k in 2:2:length(f)
                r, a = f[k], finite_number(f[k+1], ln)
                if r == objrow
                    lp.obj[v.index] += a
                elseif r in freerows
                    continue
                elseif haskey(rowindex, r)
                    c = lp.rows[rowindex[r]].coeffs
                    c[v.index] = get(c, v.index, zero(Q)) + a
                    iszero(c[v.index]) && delete!(c, v.index)
                else
                    throw(MPSParseError(ln, "unknown row \"$r\""))
                end
            end
        elseif section === :rhs || section === :ranges
            f = drop_set_name!(f, ln, section === :rhs ? "RHS" : "RANGES")
            setname = f[1]
            if section === :rhs
                rhsset === nothing && (rhsset = setname)
                setname == rhsset || continue          # only the first set is used
            else
                rngset === nothing && (rngset = setname)
                setname == rngset || continue
            end
            for k in 2:2:length(f)
                r, a = f[k], finite_number(f[k+1], ln)
                if section === :rhs && r == objrow
                    lp.objconst = -a
                elseif r in freerows || (section === :ranges && r == objrow)
                    continue
                elseif haskey(rowindex, r)
                    i = rowindex[r]
                    if section === :rhs
                        old = lp.rows[i]
                        lp.rows[i] = Row(old.name, old.coeffs, old.sense, a)
                    else
                        ranges[i] = (a, ln)
                    end
                else
                    throw(MPSParseError(ln, "unknown row \"$r\""))
                end
            end
        elseif section === :bounds
            t = uppercase(f[1])
            novalue = t in ("FR", "MI", "PL", "BV")
            t in ("UP", "LO", "FX", "FR", "MI", "PL", "BV", "LI", "UI") ||
                throw(MPSParseError(ln, t == "SC" ? "semi-continuous bounds (SC) are not supported" :
                                        "unknown bound type \"$(f[1])\""))
            # forms: TYPE SET COL [VALUE] or TYPE COL [VALUE]
            n = length(f)
            if novalue
                n in (2, 3, 4) || throw(MPSParseError(ln, "malformed $t bound"))
                setname, col = n == 2 ? ("", f[2]) : (f[2], f[3])
                n == 4 && t != "BV" && throw(MPSParseError(ln, "malformed $t bound"))
                val = nothing
            else
                n in (3, 4) || throw(MPSParseError(ln, "malformed $t bound"))
                setname, col, sval = n == 3 ? ("", f[2], f[3]) : (f[2], f[3], f[4])
                val = mps_number(sval, ln)
            end
            bndset === nothing && (bndset = setname)
            setname == bndset || isempty(setname) || continue
            haskey(colindex, col) || throw(MPSParseError(ln, "unknown column \"$col\""))
            j = colindex[col].index
            if t == "UP" || t == "UI"
                val === :neginf && throw(MPSParseError(ln, "upper bound of -inf"))
                lp.ub[j] = val === :posinf ? nothing : val
                val isa Q && val < 0 && !(j in lbset) && (lp.lb[j] = nothing)
                t == "UI" && (lp.integer[j] = true)
            elseif t == "LO" || t == "LI"
                val === :posinf && throw(MPSParseError(ln, "lower bound of +inf"))
                lp.lb[j] = val === :neginf ? nothing : val
                push!(lbset, j)
                t == "LI" && (lp.integer[j] = true)
            elseif t == "FX"
                val isa Symbol && throw(MPSParseError(ln, "fixed bound must be finite"))
                lp.lb[j] = val; lp.ub[j] = val
                push!(lbset, j)
            elseif t == "FR"
                lp.lb[j] = nothing; lp.ub[j] = nothing
                push!(lbset, j)
            elseif t == "MI"
                lp.lb[j] = nothing
                push!(lbset, j)
            elseif t == "PL"
                lp.ub[j] = nothing
            else # BV
                lp.lb[j] = zero(Q); lp.ub[j] = one(Q); lp.integer[j] = true
                push!(lbset, j)
            end
        elseif section === :name
            throw(MPSParseError(ln, "unexpected text in NAME section"))
        else
            throw(MPSParseError(ln, "data outside a section"))
        end
    end
    sawrows || throw(MPSParseError(1, "missing ROWS section"))
    isempty(objrow) && throw(MPSParseError(1, "no objective (N) row"))
    apply_ranges!(lp, ranges)
    lp
end

function set_mps_sense!(lp::LP, word::AbstractString, ln::Int)
    w = uppercase(word)
    w in ("MAX", "MAXIMIZE", "MAXIMISE") ? (lp.sense = :max) :
    w in ("MIN", "MINIMIZE", "MINIMISE") ? (lp.sense = :min) :
    throw(MPSParseError(ln, "OBJSENSE must be MAX or MIN, got \"$word\""))
end

# RHS/RANGES lines are SET ROW VALUE [ROW VALUE]; some writers omit SET.
function drop_set_name!(f, ln, what)
    length(f) in (2, 4) && return vcat([""], f)
    length(f) in (3, 5) || throw(MPSParseError(ln, "expected row/value pairs in $what"))
    f
end

function apply_ranges!(lp::LP, ranges)
    isempty(ranges) && return lp
    names = Set(r.name for r in lp.rows)
    fresh(base) = begin
        nm = base
        k = 1
        while nm in names
            k += 1
            nm = string(base, "_", k)
        end
        push!(names, nm)
        nm
    end
    rows = Row[]
    for (i, r) in enumerate(lp.rows)
        if !haskey(ranges, i)
            push!(rows, r)
            continue
        end
        R = ranges[i][1]
        if r.sense === :le
            push!(rows, r, Row(fresh(r.name * "_lo"), r.coeffs, :ge, r.rhs - abs(R)))
        elseif r.sense === :ge
            push!(rows, r, Row(fresh(r.name * "_hi"), r.coeffs, :le, r.rhs + abs(R)))
        elseif R >= 0
            push!(rows, Row(r.name, r.coeffs, :ge, r.rhs), Row(fresh(r.name * "_hi"), r.coeffs, :le, r.rhs + R))
        else
            push!(rows, Row(r.name, r.coeffs, :le, r.rhs), Row(fresh(r.name * "_lo"), r.coeffs, :ge, r.rhs + R))
        end
    end
    lp.rows = rows
    lp
end

"Read an MPS file from disk (see [`parse_mps`](@ref))."
read_mps(path::AbstractString; fixed::Bool = false) = parse_mps(read(path, String); fixed)

"""
    mps_value(q; strict = false) -> String

Render `q` for an MPS file. Rationals whose denominator has only factors 2 and
5 are written as exact decimals. Others (like 1/3) are written as `n/d`, which
[`parse_mps`](@ref) reads but other solvers do not; with `strict = true` they
throw an `ArgumentError` instead.
"""
function mps_value(q::Q; strict::Bool = false)
    d = denominator(q)
    isone(d) && return string(numerator(q))
    k2 = k5 = 0
    while iseven(d); d ÷= 2; k2 += 1; end
    while d % 5 == 0; d ÷= 5; k5 += 1; end
    if !isone(d)
        strict && throw(ArgumentError("$(fmt(q)) has no exact decimal form; MPS files hold decimals"))
        return fmt(q)
    end
    k = max(k2, k5)
    n = numerator(q) * big(10)^k ÷ denominator(q)
    s = lpad(string(abs(n)), k + 1, '0')
    string(n < 0 ? "-" : "", s[1:end-k], ".", rstrip(s[end-k+1:end], '0'))
end

"""
    write_mps(io, lp; name = "SIMPLEXACT", strict = false)

Write `lp` as free MPS that [`parse_mps`](@ref) reads back exactly (names
must not contain whitespace). See [`mps_value`](@ref) for `strict`.
"""
function write_mps(io::IO, lp::LP; name::AbstractString = "SIMPLEXACT", strict::Bool = false)
    for n in Iterators.flatten((lp.names, (r.name for r in lp.rows)))
        (isempty(n) || any(isspace, n) || startswith(n, '$')) &&
            throw(ArgumentError("name \"$n\" cannot be written to free MPS"))
    end
    nvars(lp) == 0 && throw(ArgumentError("an MPS file needs at least one column"))
    v(q) = mps_value(q; strict)
    used = Set(r.name for r in lp.rows)
    obj = "obj"
    k = 1
    while obj in used
        k += 1
        obj = "obj_$k"
    end
    println(io, "NAME          ", name)
    lp.sense === :max && println(io, "OBJSENSE\n    MAX")
    println(io, "ROWS")
    println(io, " N  ", obj)
    for r in lp.rows
        println(io, " ", r.sense === :le ? "L" : r.sense === :ge ? "G" : "E", "  ", r.name)
    end
    println(io, "COLUMNS")
    col = [Tuple{String,Q}[] for _ in 1:nvars(lp)]
    for j in 1:nvars(lp)
        iszero(lp.obj[j]) || push!(col[j], (obj, lp.obj[j]))
    end
    for r in lp.rows, j in sort!(collect(keys(r.coeffs)))
        push!(col[j], (r.name, r.coeffs[j]))
    end
    inint = false
    nmark = 0
    for j in 1:nvars(lp)
        if lp.integer[j] != inint
            inint = lp.integer[j]
            println(io, "    MARKER    'MARKER'    '", inint ? "INTORG" : "INTEND", "'")
        end
        entries = isempty(col[j]) ? [(obj, zero(Q))] : col[j]
        for (r, a) in entries
            println(io, "    ", lp.names[j], "  ", r, "  ", v(a))
        end
    end
    inint && println(io, "    MARKER    'MARKER'    'INTEND'")
    rhs = [(r.name, r.rhs) for r in lp.rows if !iszero(r.rhs)]
    iszero(lp.objconst) || pushfirst!(rhs, (obj, -lp.objconst))
    if !isempty(rhs)
        println(io, "RHS")
        for (r, a) in rhs
            println(io, "    RHS  ", r, "  ", v(a))
        end
    end
    bounds = String[]
    for j in 1:nvars(lp)
        l, u, n = lp.lb[j], lp.ub[j], lp.names[j]
        if lp.integer[j] && l == 0 && u == 1
            push!(bounds, " BV BND  $n")
        elseif l === nothing && u === nothing
            push!(bounds, " FR BND  $n")
        elseif l !== nothing && u !== nothing && l == u
            push!(bounds, " FX BND  $n  $(v(l))")
        else
            if l === nothing
                push!(bounds, " MI BND  $n")
            elseif !iszero(l) || (u !== nothing && u < 0)
                # an explicit LO keeps a negative UP from implying -inf
                push!(bounds, " LO BND  $n  $(v(l))")
            end
            u === nothing || push!(bounds, " UP BND  $n  $(v(u))")
        end
    end
    if !isempty(bounds)
        println(io, "BOUNDS")
        foreach(b -> println(io, b), bounds)
    end
    println(io, "ENDATA")
end

write_mps(lp::LP; kw...) = sprint(io -> write_mps(io, lp; kw...))

"""
    read_model(path; format = :auto) -> LP

Read an LP-format or MPS file. With `format = :auto` the extension decides
(`.mps`, `.mps.txt` and `.fixed.mps` are MPS, anything else is LP format),
falling back to sniffing for an MPS `ROWS` section.
"""
function read_model(path::AbstractString; format::Symbol = :auto)
    fmtsym = format === :auto ? guess_format(path) : format
    fmtsym === :mps ? read_mps(path) : fmtsym === :fixedmps ? read_mps(path; fixed = true) :
    fmtsym === :lp ? read_lp(path) : throw(ArgumentError("unknown format :$format"))
end

function guess_format(path::AbstractString)
    p = lowercase(path)
    endswith(p, ".fixed.mps") && return :fixedmps
    (endswith(p, ".mps") || endswith(p, ".mps.txt") || endswith(p, ".fmps")) && return :mps
    endswith(p, ".lp") && return :lp
    isfile(path) || return :lp
    for l in eachline(path)
        s = strip(l)
        (isempty(s) || startswith(s, '*') || startswith(s, '\\')) && continue
        u = uppercase(first(split(s)))
        u in ("NAME", "ROWS") && !isspace(l[1]) && return :mps
        return :lp
    end
    :lp
end
