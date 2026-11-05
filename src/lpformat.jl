# Reading and writing a CPLEX-style LP text format (a practical subset).
#
#   \ comments start with a backslash
#   maximize                      (or: max, minimize, min)
#     profit: 3 x + 2 y + 1/2
#   subject to                    (or: st, s.t., such that)
#     cap: x + y <= 4
#     x + 3 y
#       <= 6                      (a row may span lines until its sense)
#   bounds
#     x <= 3
#     -1 <= y <= 5/2
#     z free
#   general                       (integer variables; also: generals, integer, integers)
#     x
#   binary                        (0/1 variables; also: binaries, bin)
#     b
#   end

"Error raised for malformed LP text, with the 1-based line number."
struct LPParseError <: Exception
    line::Int
    msg::String
end
Base.showerror(io::IO, e::LPParseError) = print(io, "line ", e.line, ": ", e.msg)

const NUM_RE = r"^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?:/\d+)?"
const ID_RE = r"^[A-Za-z_][A-Za-z0-9_.\[\]{}#$%&~@^']*"

function tokenize(text::AbstractString, line::Int)
    toks = Tuple{Symbol,String}[]
    i = 1
    s = String(text)
    while i <= lastindex(s)
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
            continue
        end
        rest = SubString(s, i)
        if (m = match(r"^(<=|>=|=<|=>|==|<|>|=)", rest)) !== nothing
            op = m.match
            push!(toks, (:sense, op in ("<=", "=<", "<") ? "<=" : op in (">=", "=>", ">") ? ">=" : "="))
            i += ncodeunits(op)
        elseif c in ('+', '-', '*', ':')
            push!(toks, (c == ':' ? :colon : :op, string(c)))
            i += 1
        elseif (m = match(NUM_RE, rest)) !== nothing
            push!(toks, (:num, m.match))
            i += ncodeunits(m.match)
        elseif (m = match(ID_RE, rest)) !== nothing
            push!(toks, (:id, m.match))
            i += ncodeunits(m.match)
        else
            throw(LPParseError(line, "unexpected character '$c'"))
        end
    end
    toks
end

isinf_word(t) = t[1] === :id && lowercase(t[2]) in ("inf", "infinity")

"Parse a linear expression; returns (coefficients by name in order of appearance, constant)."
function parse_expr(toks, line)
    coeffs = Pair{String,Q}[]
    con = zero(Q)
    i = 1
    n = length(toks)
    n == 0 && return coeffs, con
    while i <= n
        sign = one(Q)
        sawsign = false
        while i <= n && toks[i][1] === :op && toks[i][2] in ("+", "-")
            toks[i][2] == "-" && (sign = -sign)
            sawsign = true
            i += 1
        end
        (i == 1 || sawsign) || throw(LPParseError(line, "expected '+' or '-' before \"$(toks[i][2])\""))
        i <= n || throw(LPParseError(line, "expression ends with a sign"))
        coef = nothing
        if toks[i][1] === :num
            coef = parse_q(toks[i][2])
            i += 1
            if i <= n && toks[i] == (:op, "*")
                i += 1
                (i <= n && toks[i][1] === :id) || throw(LPParseError(line, "expected a variable after '*'"))
            end
        end
        if i <= n && toks[i][1] === :id
            isinf_word(toks[i]) && throw(LPParseError(line, "infinity is only allowed in bounds"))
            push!(coeffs, toks[i][2] => sign * (coef === nothing ? one(Q) : coef))
            i += 1
        elseif coef !== nothing
            con += sign * coef
        else
            throw(LPParseError(line, "unexpected \"$(toks[i][2])\""))
        end
    end
    coeffs, con
end

function strip_name(toks)
    if length(toks) >= 2 && toks[1][1] === :id && toks[2][1] === :colon
        return toks[1][2], toks[3:end]
    end
    any(t -> t[1] === :colon, toks) && return nothing, nothing
    "", toks
end

const SECTION_WORDS = Dict(
    "maximize" => :max, "maximise" => :max, "maximum" => :max, "max" => :max,
    "minimize" => :min, "minimise" => :min, "minimum" => :min, "min" => :min,
    "subject to" => :st, "such that" => :st, "st" => :st, "s.t." => :st,
    "bounds" => :bounds, "bound" => :bounds,
    "general" => :general, "generals" => :general, "gen" => :general,
    "integer" => :general, "integers" => :general,
    "binary" => :binary, "binaries" => :binary, "bin" => :binary,
    "end" => :end)

"""
    parse_lp(text) -> LP

Parse LP-format text (see the top of `lpformat.jl` for the grammar). Variables
default to `0 <= x < inf`, as in the CPLEX format.
"""
function parse_lp(text::AbstractString)
    lp = LP(:min)
    index = Dict{String,Var}()
    var!(name) = get!(() -> add_variable!(lp, name), index, name)
    section = :none
    sawobj = false
    obj_toks = Tuple{Symbol,String}[]
    obj_line = 0
    pending = Tuple{Symbol,String}[]
    pending_line = 0
    ended = false

    function flush_objective()
        isempty(obj_toks) && return
        name, toks = strip_name(obj_toks)
        toks === nothing && throw(LPParseError(obj_line, "malformed objective name"))
        any(t -> t[1] === :sense, toks) && throw(LPParseError(obj_line, "the objective cannot contain a comparison"))
        coeffs, con = parse_expr(toks, obj_line)
        set_objective!(lp, [var!(k) => v for (k, v) in coeffs]; constant = con)
        empty!(obj_toks)
    end
    function add_row(toks, line)
        name, toks = strip_name(toks)
        toks === nothing && throw(LPParseError(line, "malformed row name"))
        k = findall(t -> t[1] === :sense, toks)
        length(k) == 1 || throw(LPParseError(line, length(k) == 0 ? "row has no comparison" :
                                                  "ranged rows are not supported; write two rows"))
        k = k[1]
        lc, lk = parse_expr(toks[1:k-1], line)
        rc, rk = parse_expr(toks[k+1:end], line)
        isempty(lc) && isempty(rc) && throw(LPParseError(line, "row has no variables"))
        terms = vcat([var!(n) => v for (n, v) in lc], [var!(n) => -v for (n, v) in rc])
        sense = toks[k][2] == "<=" ? :le : toks[k][2] == ">=" ? :ge : :eq
        name in (r.name for r in lp.rows) && throw(LPParseError(line, "duplicate row name \"$name\""))
        add_constraint!(lp, terms, sense, rk - lk; name)
    end
    function add_bound(toks, line)
        bad() = throw(LPParseError(line, "malformed bound"))
        num(t) = isinf_word(t) ? nothing : t[1] === :num ? parse_q(t[2]) : bad()
        # fold a leading sign into the following number or infinity
        vals = Tuple{Symbol,Any}[]
        i = 1
        while i <= length(toks)
            t = toks[i]
            if t[1] === :op && t[2] in ("+", "-") && i < length(toks) && (toks[i+1][1] === :num || isinf_word(toks[i+1]))
                v = num(toks[i+1])
                push!(vals, (:val, v === nothing ? (t[2] == "-" ? :neginf : :posinf) : (t[2] == "-" ? -v : v)))
                i += 2
            elseif t[1] === :num || isinf_word(t)
                v = num(t)
                push!(vals, (:val, v === nothing ? :posinf : v))
                i += 1
            elseif t[1] === :id
                push!(vals, (:var, t[2]))
                i += 1
            elseif t[1] === :sense
                push!(vals, (:sense, t[2]))
                i += 1
            else
                bad()
            end
        end
        setlb(v, x) = (lp.lb[v.index] = x === :neginf ? nothing : x === :posinf ? bad() : x)
        setub(v, x) = (lp.ub[v.index] = x === :posinf ? nothing : x === :neginf ? bad() : x)
        kinds = [v[1] for v in vals]
        if kinds == [:var, :var] && lowercase(vals[2][2]) == "free"
            v = var!(vals[1][2]); lp.lb[v.index] = nothing; lp.ub[v.index] = nothing
        elseif kinds == [:var, :sense, :val]
            v = var!(vals[1][2]); op = vals[2][2]; x = vals[3][2]
            op == "<=" ? setub(v, x) : op == ">=" ? setlb(v, x) : (setlb(v, x); setub(v, x))
        elseif kinds == [:val, :sense, :var]
            v = var!(vals[3][2]); op = vals[2][2]; x = vals[1][2]
            op == "<=" ? setlb(v, x) : op == ">=" ? setub(v, x) : (setlb(v, x); setub(v, x))
        elseif kinds == [:val, :sense, :var, :sense, :val] && vals[2][2] == vals[4][2] && vals[2][2] != "="
            v = var!(vals[3][2])
            lo, hi = vals[2][2] == "<=" ? (vals[1][2], vals[5][2]) : (vals[5][2], vals[1][2])
            setlb(v, lo); setub(v, hi)
        else
            bad()
        end
    end

    function parse_line!(line, ln)
        toks = tokenize(line, ln)
        if section in (:max, :min)
            isempty(obj_toks) && (obj_line = ln)
            append!(obj_toks, toks)
        elseif section === :st
            isempty(pending) && (pending_line = ln)
            append!(pending, toks)
            k = findlast(t -> t[1] === :sense, pending)
            if k !== nothing && k < length(pending)
                add_row(pending, pending_line)
                empty!(pending)
            end
        elseif section === :bounds
            add_bound(toks, ln)
        elseif section in (:general, :binary)
            for t in toks
                t[1] === :id || throw(LPParseError(ln, "expected variable names"))
                v = var!(t[2])
                lp.integer[v.index] = true
                if section === :binary
                    lp.lb[v.index] = zero(Q)
                    lp.ub[v.index] = one(Q)
                end
            end
        else
            throw(LPParseError(ln, "expected a section keyword such as 'maximize' or 'minimize'"))
        end
    end

    for (ln, raw) in enumerate(split(text, '\n'))
        line = strip(first(split(raw, '\\'; limit = 2)))
        isempty(line) && continue
        ended && throw(LPParseError(ln, "text after 'end'"))
        key = lowercase(replace(line, r"\s+" => " "))
        if haskey(SECTION_WORDS, key)
            sec = SECTION_WORDS[key]
            isempty(pending) || throw(LPParseError(pending_line, "row has no comparison"))
            if sec in (:max, :min)
                sawobj && throw(LPParseError(ln, "more than one objective section"))
                sawobj = true
                lp.sense = sec
                obj_line = ln
            elseif section in (:max, :min)
                flush_objective()
            end
            sec === :end && (ended = true)
            section = sec
            continue
        end
        try
            parse_line!(line, ln)
        catch e
            e isa ArgumentError || rethrow()
            throw(LPParseError(ln, e.msg))
        end
    end
    isempty(pending) || throw(LPParseError(pending_line, "row has no comparison"))
    section in (:max, :min) && flush_objective()
    sawobj || throw(LPParseError(1, "missing 'maximize' or 'minimize' section"))
    lp
end


"Read an LP file from disk."
read_lp(path::AbstractString) = parse_lp(read(path, String))

function write_expr(io::IO, lp::LP, terms)
    first = true
    for (j, a) in terms
        iszero(a) && continue
        if first
            print(io, a < 0 ? "- " : "")
        else
            print(io, a < 0 ? " - " : " + ")
        end
        isone(abs(a)) || print(io, fmt(abs(a)), " ")
        print(io, lp.names[j])
        first = false
    end
    first
end

"""
    write_lp(io, lp)

Write `lp` in the LP format that [`parse_lp`](@ref) reads back exactly.
"""
function write_lp(io::IO, lp::LP)
    println(io, lp.sense === :max ? "maximize" : "minimize")
    print(io, "  obj: ")
    empty = write_expr(io, lp, ((j, lp.obj[j]) for j in 1:nvars(lp)))
    if !iszero(lp.objconst)
        print(io, empty ? (lp.objconst < 0 ? "- " : "") : (lp.objconst < 0 ? " - " : " + "), fmt(abs(lp.objconst)))
    elseif empty
        print(io, "0")
    end
    println(io)
    println(io, "subject to")
    for r in lp.rows
        print(io, "  ", r.name, ": ")
        if write_expr(io, lp, sort!(collect(r.coeffs)))
            # an all-zero row still needs a variable for the parser
            print(io, "0 ", lp.names[1])
        end
        println(io, " ", r.sense === :le ? "<=" : r.sense === :ge ? ">=" : "=", " ", fmt(r.rhs))
    end
    println(io, "bounds")
    for j in 1:nvars(lp)
        l, u = lp.lb[j], lp.ub[j]
        n = lp.names[j]
        if l === nothing && u === nothing
            println(io, "  ", n, " free")
        elseif l === nothing
            println(io, "  -inf <= ", n, " <= ", fmt(u))
        elseif u === nothing
            iszero(l) ? println(io, "  ", n, " >= 0") : println(io, "  ", n, " >= ", fmt(l))
        else
            println(io, "  ", fmt(l), " <= ", n, " <= ", fmt(u))
        end
    end
    ints = [lp.names[j] for j in 1:nvars(lp) if lp.integer[j]]
    if !isempty(ints)
        println(io, "general")
        println(io, "  ", join(ints, " "))
    end
    println(io, "end")
end

write_lp(lp::LP) = sprint(write_lp, lp)
