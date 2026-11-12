# Feasibility repair: when an LP has no feasible point, find the point that
# violates its rows (and optionally its bounds) by the least weighted total.

"""
Result of [`feasibility_repair`](@ref).

* `status`: `:optimal` when a least-violation point was found, `:infeasible`
  when even the relaxed problem has no point (contradictory bounds that may
  not be relaxed, or integrality that the bounds rule out), `:unbounded` when
  `optimize = true` and the original objective is unbounded among the
  least-violation points, `:node_limit` when branch and bound stopped early;
* `violation`: the least weighted total violation (zero for a feasible LP);
* `x`: the repaired point;
* `row_violation[i]`: how far `x` misses row `i` (zero when it holds);
* `lower_violation[j]`, `upper_violation[j]`: how far `x[j]` is below its lower
  or above its upper bound (always zero unless bounds were relaxed);
* `objective`: the original objective at `x`;
* `model`, `solution`: the elastic LP that was solved last and its solution,
  whose certificate proves `violation` is the least possible.
"""
struct Repair
    status::Symbol
    violation::Union{Nothing,Q}
    x::Vector{Q}
    row_violation::Vector{Q}
    lower_violation::Vector{Q}
    upper_violation::Vector{Q}
    objective::Union{Nothing,Q}
    model::LP
    solution::Union{Solution,MIPSolution}
    weights::Vector{Q}
    bound_weight::Union{Nothing,Q}
    optimize::Bool
end

"Amount by which `x` violates row `r`."
function row_violation(r::Row, x::AbstractVector)
    v = row_activity(r, x)
    r.sense === :le ? max(zero(Q), v - r.rhs) : r.sense === :ge ? max(zero(Q), r.rhs - v) : abs(v - r.rhs)
end

lower_violation(lp::LP, j, x) = lp.lb[j] === nothing ? zero(Q) : max(zero(Q), lp.lb[j] - x[j])
upper_violation(lp::LP, j, x) = lp.ub[j] === nothing ? zero(Q) : max(zero(Q), x[j] - lp.ub[j])

"""
The elastic version of `lp`: every row gets non-negative slack variables
(`s` for a `<=` or `>=` row, `p` and `n` for an equality) whose weighted sum is
minimized. With `bound_weight` set, finite bounds become elastic rows too and
the variables keep only their integrality.
"""
function elastic_model(lp::LP, weights::Vector{Q}, bound_weight::Union{Nothing,Q})
    m = LP(:min)
    for j in 1:nvars(lp)
        if bound_weight === nothing
            add_variable!(m, lp.names[j]; lb = lp.lb[j], ub = lp.ub[j], integer = lp.integer[j])
        else
            add_variable!(m, lp.names[j]; lb = nothing, integer = lp.integer[j])
        end
    end
    taken = Set(m.names)
    fresh(base) = begin
        nm = base; k = 1
        while nm in taken
            k += 1; nm = string(base, "_", k)
        end
        push!(taken, nm); nm
    end
    obj = Pair{Var,Q}[]
    rowname(i) = lp.rows[i].name
    for (i, r) in enumerate(lp.rows)
        terms = Pair{Var,Q}[Var(j) => a for (j, a) in r.coeffs]
        w = weights[i]
        if r.sense === :le
            s = add_variable!(m, fresh("viol_" * rowname(i)))
            push!(terms, s => -1); push!(obj, s => w)
        elseif r.sense === :ge
            s = add_variable!(m, fresh("viol_" * rowname(i)))
            push!(terms, s => 1); push!(obj, s => w)
        else
            p = add_variable!(m, fresh("under_" * rowname(i)))
            n = add_variable!(m, fresh("over_" * rowname(i)))
            push!(terms, p => 1, n => -1); push!(obj, p => w, n => w)
        end
        add_constraint!(m, terms, r.sense, r.rhs; name = r.name)
    end
    if bound_weight !== nothing
        rownames = Set(r.name for r in m.rows)
        rfresh(base) = begin
            nm = base; k = 1
            while nm in rownames
                k += 1; nm = string(base, "_", k)
            end
            push!(rownames, nm); nm
        end
        for j in 1:nvars(lp)
            if lp.lb[j] !== nothing
                s = add_variable!(m, fresh("below_" * lp.names[j]))
                push!(obj, s => bound_weight)
                add_constraint!(m, [Var(j) => 1, s => 1], :ge, lp.lb[j]; name = rfresh("lb_" * lp.names[j]))
            end
            if lp.ub[j] !== nothing
                s = add_variable!(m, fresh("above_" * lp.names[j]))
                push!(obj, s => bound_weight)
                add_constraint!(m, [Var(j) => 1, s => -1], :le, lp.ub[j]; name = rfresh("ub_" * lp.names[j]))
            end
        end
    end
    set_objective!(m, obj)
    m
end

"""
    feasibility_repair(lp; weights = nothing, relax_bounds = false, bound_weight = 1,
                       optimize = false, rule = :auto, node_limit = 100_000) -> Repair

Find a point that violates the rows of `lp` by the least weighted total, the
classic elastic (or "phase 1 with penalties") repair of an infeasible model.
`weights` gives one positive weight per row, or a `Dict` from row name to
weight (missing rows get 1). With `relax_bounds = true` variable bounds may be
violated too, at `bound_weight` per unit. Integer variables stay integer (the
elastic model is then solved by branch and bound).

With `optimize = true` a second solve optimizes the original objective among
all least-violation points. A feasible LP repairs with violation zero, and
with `optimize` the repair is then its optimum.

Every answer comes with exact certificates; [`verify_repair`](@ref) checks
them against `lp` without trusting the repair code.
"""
function feasibility_repair(lp::LP; weights = nothing, relax_bounds::Bool = false, bound_weight = 1,
        optimize::Bool = false, rule::Symbol = :auto, node_limit::Int = 100_000)
    w = repair_weights(lp, weights)
    bw = relax_bounds ? tonum(bound_weight) : nothing
    bw !== nothing && bw <= 0 && throw(ArgumentError("bound_weight must be positive"))
    m = elastic_model(lp, w, bw)
    sol, vmin = solve_model(m; rule, node_limit)
    if sol.status !== :optimal
        return Repair(sol.status, nothing, Q[], Q[], Q[], Q[], nothing, m, sol, w, bw, optimize)
    end
    if optimize
        m2 = elastic_model(lp, w, bw)
        add_constraint!(m2, [Var(j) => a for (j, a) in enumerate(m.obj) if !iszero(a)], :le, vmin;
            name = repair_rowname(m2))
        m2.sense = lp.sense
        fill!(m2.obj, zero(Q))
        m2.obj[1:nvars(lp)] .= lp.obj
        m2.objconst = lp.objconst
        sol2, _ = solve_model(m2; rule, node_limit)
        if sol2.status !== :optimal
            return Repair(sol2.status, vmin, Q[], Q[], Q[], Q[], nothing, m2, sol2, w, bw, true)
        end
        m, sol = m2, sol2
    end
    x = sol.x[1:nvars(lp)]
    rv = [row_violation(r, x) for r in lp.rows]
    lv = [lower_violation(lp, j, x) for j in 1:nvars(lp)]
    uv = [upper_violation(lp, j, x) for j in 1:nvars(lp)]
    Repair(:optimal, vmin, x, rv, lv, uv, objective_value(lp, x), m, sol, w, bw, optimize)
end

function repair_rowname(m::LP)
    names = Set(r.name for r in m.rows)
    nm = "least_violation"; k = 1
    while nm in names
        k += 1; nm = "least_violation_$k"
    end
    nm
end

function solve_model(m::LP; rule, node_limit)
    if any(m.integer)
        s = solve_mip(m; rule, node_limit)
        s, s.objective
    else
        s = solve(m; rule)
        s, s.objective
    end
end

function repair_weights(lp::LP, weights)
    w = if weights === nothing
        fill(one(Q), nrows(lp))
    elseif weights isa AbstractDict
        names = Set(r.name for r in lp.rows)
        for k in keys(weights)
            String(k) in names || throw(ArgumentError("no row named \"$k\""))
        end
        [tonum(get(weights, r.name, 1)) for r in lp.rows]
    else
        length(weights) == nrows(lp) || throw(ArgumentError("need one weight per row"))
        tonum.(collect(weights))
    end
    all(>(0), w) || throw(ArgumentError("weights must be positive"))
    Vector{Q}(w)
end

"""
    check_repair(lp, r) -> Union{Nothing,String}

Check a [`Repair`](@ref) against `lp`: the elastic model is rebuilt from `lp`
(so a tampered model is caught), its solution's certificate must verify, the
reported violations must be what `x` actually does to `lp`, and their weighted
total must be the proven least violation.
"""
function check_repair(lp::LP, r::Repair)
    r.status === :optimal || return "status is $(r.status), nothing to check"
    m = elastic_model(lp, r.weights, r.bound_weight)
    if r.optimize
        add_constraint!(m, [Var(j) => a for (j, a) in enumerate(m.obj) if !iszero(a)], :le, r.violation;
            name = repair_rowname(m))
        m.sense = lp.sense
        fill!(m.obj, zero(Q))
        m.obj[1:nvars(lp)] .= lp.obj
        m.objconst = lp.objconst
    end
    write_lp(m) == write_lp(r.model) || return "the elastic model does not match the problem"
    reason = r.solution isa MIPSolution ? check_mip_certificate(m, r.solution) : check_certificate(m, r.solution)
    reason === nothing || return "elastic model certificate: " * reason
    r.x == r.solution.x[1:nvars(lp)] || return "x does not match the elastic solution"
    for j in 1:nvars(lp)
        lp.integer[j] && !isone(denominator(r.x[j])) && return "integer variable $j is fractional"
    end
    r.row_violation == [row_violation(row, r.x) for row in lp.rows] || return "row violations are misreported"
    r.lower_violation == [lower_violation(lp, j, r.x) for j in 1:nvars(lp)] || return "lower-bound violations are misreported"
    r.upper_violation == [upper_violation(lp, j, r.x) for j in 1:nvars(lp)] || return "upper-bound violations are misreported"
    if r.bound_weight === nothing
        (all(iszero, r.lower_violation) && all(iszero, r.upper_violation)) || return "bounds violated but not relaxable"
    end
    total = sum(r.weights .* r.row_violation; init = zero(Q)) +
        (r.bound_weight === nothing ? zero(Q) : r.bound_weight * (sum(r.lower_violation; init = zero(Q)) + sum(r.upper_violation; init = zero(Q))))
    if r.optimize
        total <= r.violation || return "x violates more than the least violation"
        # the least violation itself is proven by re-solving the elastic model
        m1 = elastic_model(lp, r.weights, r.bound_weight)
        s1, v1 = solve_model(m1; rule = :bland, node_limit = 100_000)
        s1.status === :optimal && v1 == r.violation || return "the least violation is not $(fmt(r.violation))"
        c1 = s1 isa MIPSolution ? check_mip_certificate(m1, s1) : check_certificate(m1, s1)
        c1 === nothing || return "least-violation certificate: " * c1
    else
        total == r.violation || return "weighted violation $(fmt(total)) is not the proven least $(fmt(r.violation))"
    end
    r.objective == objective_value(lp, r.x) || return "objective is misreported"
    nothing
end

"`true` when `r` is a verified least-violation repair of `lp`."
verify_repair(lp::LP, r::Repair) = check_repair(lp, r) === nothing
