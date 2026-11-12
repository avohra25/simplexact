# Irreducible infeasible subsystems: a smallest-possible explanation of why an
# LP has no feasible point.

"""
An irreducible infeasible subsystem (IIS) of an [`LP`](@ref): a set of rows
and variable bounds that has no feasible point on its own, while dropping any
single member makes it feasible.

* `rows`: indices of the member rows;
* `lower`, `upper`: variables whose lower or upper bound is a member (all
  other bounds are dropped, so those variables are free in the subsystem);
* `farkas`: one multiplier per member row, a certificate that the subsystem is
  infeasible (empty when the IIS is a variable with `lb > ub`);
* `witnesses`: for every member, a point that satisfies all other members,
  proving that the member cannot be dropped. Keys are `(:row, i)`,
  `(:lower, j)` or `(:upper, j)`.

Integrality is ignored: this explains infeasibility of the LP relaxation.
"""
struct IIS
    rows::Vector{Int}
    lower::Vector{Int}
    upper::Vector{Int}
    farkas::Vector{Q}
    witnesses::Dict{Tuple{Symbol,Int},Vector{Q}}
end

nmembers(iis::IIS) = length(iis.rows) + length(iis.lower) + length(iis.upper)

"The LP keeping only rows `rows` and the bounds in `lower`/`upper`, with a zero objective."
function subsystem(lp::LP, rows, lower, upper)
    sub = LP(:min)
    for j in 1:nvars(lp)
        add_variable!(sub, lp.names[j]; lb = j in lower ? lp.lb[j] : nothing,
            ub = j in upper ? lp.ub[j] : nothing)
    end
    for i in rows
        r = lp.rows[i]
        push!(sub.rows, Row(r.name, copy(r.coeffs), r.sense, r.rhs))
    end
    sub
end

# The rows and bounds a Farkas certificate actually uses: rows with a non-zero
# multiplier, and the bound each variable of the combined row is pushed against.
function farkas_support(lp::LP, rows, u)
    keep = [rows[k] for k in eachindex(rows) if !iszero(u[k])]
    w = zeros(Q, nvars(lp))
    for (k, i) in enumerate(rows), (j, a) in lp.rows[i].coeffs
        w[j] += u[k] * a
    end
    keep, Set(j for j in 1:nvars(lp) if w[j] > 0), Set(j for j in 1:nvars(lp) if w[j] < 0)
end

"""
    find_iis(lp; rule = :auto) -> IIS

Find an irreducible infeasible subsystem of an infeasible LP (integrality is
ignored). Throws an `ArgumentError` when the LP relaxation is feasible.

The search starts from the rows and bounds the solver's Farkas certificate
uses, then runs a deletion filter: each member is dropped in turn, and stays
dropped if what is left is still infeasible (shrinking the candidate set to
the new certificate's support). Members that cannot be dropped get the
feasible point that proves it. Every step is an exact solve, so the result is
irreducible, not just small. Check it with [`verify_iis`](@ref).
"""
function find_iis(lp::LP; rule::Symbol = :auto)
    j = conflicting_bounds(lp)
    if j != 0
        w = Dict{Tuple{Symbol,Int},Vector{Q}}()
        x = zeros(Q, nvars(lp))
        x[j] = lp.lb[j]; w[(:upper, j)] = copy(x)
        x[j] = lp.ub[j]; w[(:lower, j)] = copy(x)
        return IIS(Int[], [j], [j], Q[], w)
    end
    sol = solve(lp; rule)
    sol.status === :infeasible || throw(ArgumentError("the LP relaxation is feasible ($(sol.status)); there is no IIS"))
    rows, lower, upper = farkas_support(lp, collect(1:nrows(lp)), sol.farkas)
    witnesses = Dict{Tuple{Symbol,Int},Vector{Q}}()
    candidates = vcat([(:row, i) for i in rows], [(:lower, j) for j in sort!(collect(lower))],
        [(:upper, j) for j in sort!(collect(upper))])
    for m in candidates
        member(m, rows, lower, upper) || continue          # already filtered out
        r2, l2, u2 = without(m, rows, lower, upper)
        s = solve(subsystem(lp, r2, l2, u2); rule)
        if s.status === :infeasible
            rows, lower, upper = farkas_support(lp, r2, s.farkas)
        else
            witnesses[m] = s.x
        end
    end
    # all members survived the filter; one last solve gives the certificate
    s = solve(subsystem(lp, rows, lower, upper); rule)
    s.status === :infeasible || error("internal error: IIS candidate became feasible")
    IIS(rows, sort!(collect(lower)), sort!(collect(upper)), s.farkas,
        Dict(m => witnesses[m] for m in keys(witnesses) if member(m, rows, lower, upper)))
end

member(m, rows, lower, upper) = m[1] === :row ? m[2] in rows : m[1] === :lower ? m[2] in lower : m[2] in upper

function without(m, rows, lower, upper)
    m[1] === :row && return filter(!=(m[2]), rows), lower, upper
    m[1] === :lower && return rows, setdiff(lower, (m[2],)), upper
    rows, lower, setdiff(upper, (m[2],))
end

"""
    check_iis(lp, iis) -> Union{Nothing,String}

Independently check that `iis` is an irreducible infeasible subsystem of
`lp`: its Farkas certificate proves the subsystem infeasible, and each
member's witness point satisfies every other member. Returns `nothing` or a
reason.
"""
function check_iis(lp::LP, iis::IIS)
    nv = nvars(lp)
    all(i -> 1 <= i <= nrows(lp), iis.rows) || return "row index out of range"
    allunique(iis.rows) && allunique(iis.lower) && allunique(iis.upper) || return "repeated member"
    for j in iis.lower
        1 <= j <= nv && lp.lb[j] !== nothing || return "variable $j has no finite lower bound"
    end
    for j in iis.upper
        1 <= j <= nv && lp.ub[j] !== nothing || return "variable $j has no finite upper bound"
    end
    nmembers(iis) > 0 || return "empty subsystem"
    # infeasibility
    if isempty(iis.rows) && length(iis.lower) == 1 && iis.lower == iis.upper &&
            lp.lb[iis.lower[1]] > lp.ub[iis.lower[1]]
        # a single contradictory pair of bounds needs no multipliers
    else
        sub = subsystem(lp, iis.rows, Set(iis.lower), Set(iis.upper))
        conflicting_bounds(sub) == 0 || return "the subsystem has contradictory bounds but is not just that pair"
        reason = check_certificate(sub, Solution(:infeasible; farkas = iis.farkas))
        reason === nothing || return "subsystem is not proven infeasible: " * reason
    end
    # irreducibility
    members = vcat([(:row, i) for i in iis.rows], [(:lower, j) for j in iis.lower], [(:upper, j) for j in iis.upper])
    for m in members
        x = get(iis.witnesses, m, nothing)
        x === nothing && return "no witness for $(m[1]) $(m[2])"
        length(x) == nv || return "witness for $(m[1]) $(m[2]) has the wrong length"
        for i in iis.rows
            (m == (:row, i) || row_satisfied(lp.rows[i], x)) || return "witness for $(m[1]) $(m[2]) violates row $i"
        end
        for j in iis.lower
            (m == (:lower, j) || x[j] >= lp.lb[j]) || return "witness for $(m[1]) $(m[2]) violates the lower bound of variable $j"
        end
        for j in iis.upper
            (m == (:upper, j) || x[j] <= lp.ub[j]) || return "witness for $(m[1]) $(m[2]) violates the upper bound of variable $j"
        end
    end
    nothing
end

"`true` when `iis` is a verified irreducible infeasible subsystem of `lp`."
verify_iis(lp::LP, iis::IIS) = check_iis(lp, iis) === nothing
