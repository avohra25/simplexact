# Problem model: variables with exact rational bounds, linear rows, objective.

"Exact rational number type used everywhere in the solver."
const Q = Rational{BigInt}

"A bound is a rational number, or `nothing` for 'no bound'."
const Bound = Union{Nothing,Q}

"""
    parse_q(s)

Parse an exact rational from text: `"3"`, `"-7/2"`, `"0.125"`, `"1.5e-3"`.
"""
function parse_q(s::AbstractString)
    s = strip(s)
    if occursin('/', s)
        a, b = split(s, '/'; limit = 2)
        d = parse_q(b)
        iszero(d) && throw(ArgumentError("zero denominator in \"$s\""))
        return parse_q(a) / d
    end
    m = match(r"^([+-]?)(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$", s)
    m === nothing && throw(ArgumentError("invalid number: \"$s\""))
    ip = m[2]
    fp = m[3] === nothing ? "" : String(m[3])
    isempty(ip) && isempty(fp) && throw(ArgumentError("invalid number: \"$s\""))
    digits = ip * fp
    num = parse(BigInt, isempty(digits) ? "0" : digits)
    q = Q(num, big(10)^length(fp))
    if m[4] !== nothing
        e = parse(Int, m[4])
        q *= e >= 0 ? Q(big(10)^e) : Q(1, big(10)^(-e))
    end
    m[1] == "-" ? -q : q
end

"Convert integers, rationals, decimal strings and floats (via their shortest decimal form) to `Q`."
tonum(x::Q) = x
tonum(x::Integer) = Q(x)
tonum(x::Rational) = Q(x.num, x.den)
tonum(x::AbstractString) = parse_q(x)
function tonum(x::AbstractFloat)
    isfinite(x) || throw(ArgumentError("non-finite coefficient $x; use `nothing` for an infinite bound"))
    parse_q(repr(x))
end

tobound(::Nothing) = nothing
tobound(x) = tonum(x)

"Format a rational as `n` or `n/d`."
fmt(q::Rational) = isone(q.den) ? string(q.num) : string(q.num, "/", q.den)

"Approximate decimal rendering of a rational with `digits` places."
function decimal(q::Rational; digits::Int = 6)
    neg = q < 0
    n = abs(numerator(q)) * big(10)^digits
    r = div(n + div(denominator(q), 2), denominator(q))   # round half up
    s = string(r)
    s = lpad(s, digits + 1, '0')
    ip, fp = s[1:end-digits], s[end-digits+1:end]
    string(neg && r != 0 ? "-" : "", ip, digits > 0 ? "." * fp : "")
end

"A handle to a variable of an [`LP`](@ref)."
struct Var
    index::Int
end

"A linear row `sum(coeffs[j] * x[j]) <sense> rhs` with `sense` one of `:le`, `:ge`, `:eq`."
struct Row
    name::String
    coeffs::Dict{Int,Q}
    sense::Symbol
    rhs::Q
end

"""
    LP(sense = :min)

A linear (or mixed-integer) program with exact rational data. `sense` is
`:min` or `:max`. Build it with [`add_variable!`](@ref),
[`set_objective!`](@ref) and [`add_constraint!`](@ref), then call [`solve`](@ref).
"""
mutable struct LP
    sense::Symbol
    names::Vector{String}
    lb::Vector{Bound}
    ub::Vector{Bound}
    integer::Vector{Bool}
    obj::Vector{Q}
    objconst::Q
    rows::Vector{Row}
end

function LP(sense::Symbol = :min)
    sense in (:min, :max) || throw(ArgumentError("sense must be :min or :max, got :$sense"))
    LP(sense, String[], Bound[], Bound[], Bool[], Q[], zero(Q), Row[])
end

nvars(lp::LP) = length(lp.names)
nrows(lp::LP) = length(lp.rows)
Base.getindex(lp::LP, name::AbstractString) = begin
    i = findfirst(==(String(name)), lp.names)
    i === nothing && throw(KeyError(name))
    Var(i)
end

"""
    add_variable!(lp, name; lb = 0, ub = nothing, integer = false) -> Var

Add a variable. `lb = nothing` means no lower bound (a free variable when `ub`
is also `nothing`).
"""
function add_variable!(lp::LP, name::AbstractString; lb = 0, ub = nothing, integer::Bool = false)
    any(==(String(name)), lp.names) && throw(ArgumentError("duplicate variable name \"$name\""))
    push!(lp.names, String(name))
    push!(lp.lb, tobound(lb))
    push!(lp.ub, tobound(ub))
    push!(lp.integer, integer)
    push!(lp.obj, zero(Q))
    Var(nvars(lp))
end

function collect_terms(lp::LP, terms)
    d = Dict{Int,Q}()
    for (v, a) in terms
        v isa Var || throw(ArgumentError("expected a Var, got $(typeof(v))"))
        1 <= v.index <= nvars(lp) || throw(ArgumentError("variable $(v.index) does not belong to this LP"))
        d[v.index] = get(d, v.index, zero(Q)) + tonum(a)
    end
    filter!(p -> !iszero(p.second), d)
    d
end

"""
    set_objective!(lp, terms; sense = lp.sense, constant = 0)

Set the objective to `sum(coef * var for (var, coef) in terms) + constant`.
"""
function set_objective!(lp::LP, terms; sense::Symbol = lp.sense, constant = 0)
    sense in (:min, :max) || throw(ArgumentError("sense must be :min or :max"))
    lp.sense = sense
    fill!(lp.obj, zero(Q))
    for (j, a) in collect_terms(lp, terms)
        lp.obj[j] = a
    end
    lp.objconst = tonum(constant)
    lp
end

const SENSES = Dict{Any,Symbol}(:le => :le, :ge => :ge, :eq => :eq,
    "<=" => :le, "<" => :le, ">=" => :ge, ">" => :ge, "=" => :eq, "==" => :eq,
    :(<=) => :le, :(>=) => :ge, :(==) => :eq)

"""
    add_constraint!(lp, terms, sense, rhs; name = "") -> Int

Add the row `sum(coef * var) sense rhs`, where `sense` is `:le`, `:ge`, `:eq`
(or the strings `"<="`, `">="`, `"="`). Returns the row index.
"""
function add_constraint!(lp::LP, terms, sense, rhs; name::AbstractString = "")
    s = get(SENSES, sense, nothing)
    s === nothing && throw(ArgumentError("unknown constraint sense $(repr(sense))"))
    nm = isempty(name) ? "c$(nrows(lp) + 1)" : String(name)
    push!(lp.rows, Row(nm, collect_terms(lp, terms), s, tonum(rhs)))
    nrows(lp)
end

"Objective value of `x` (a vector of rationals, one per variable)."
objective_value(lp::LP, x::AbstractVector) =
    lp.objconst + sum((lp.obj[j] * x[j] for j in 1:nvars(lp)); init = zero(Q))

"Left-hand side of row `r` evaluated at `x`."
row_activity(r::Row, x::AbstractVector) =
    sum((a * x[j] for (j, a) in r.coeffs); init = zero(Q))

function row_satisfied(r::Row, x::AbstractVector)
    v = row_activity(r, x)
    r.sense === :le ? v <= r.rhs : r.sense === :ge ? v >= r.rhs : v == r.rhs
end

"Is `x` feasible for every bound and row (integrality is checked separately)?"
function is_feasible(lp::LP, x::AbstractVector)
    length(x) == nvars(lp) || return false
    for j in 1:nvars(lp)
        lp.lb[j] !== nothing && x[j] < lp.lb[j] && return false
        lp.ub[j] !== nothing && x[j] > lp.ub[j] && return false
    end
    all(r -> row_satisfied(r, x), lp.rows)
end
