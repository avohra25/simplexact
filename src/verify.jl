# Independent exact verification of a Solution against the original LP.
# Nothing here touches the simplex tableau: it only uses the rows, bounds and
# objective, so a bug in the solver cannot make a wrong answer verify.

"A column-wise view of the rows: `col[j] = [(row index, coefficient), ...]`."
function columns(lp::LP)
    col = [Tuple{Int,Q}[] for _ in 1:nvars(lp)]
    for (i, r) in enumerate(lp.rows)
        for (j, a) in r.coeffs
            push!(col[j], (i, a))
        end
    end
    col
end

"Row-combination `sum_i u[i] * A[i, j]` for every variable `j`."
function combine(lp::LP, u::AbstractVector)
    w = zeros(Q, nvars(lp))
    for (i, r) in enumerate(lp.rows), (j, a) in r.coeffs
        w[j] += u[i] * a
    end
    w
end

"""
    check_certificate(lp, sol) -> Union{Nothing,String}

Return `nothing` if `sol` carries a valid exact certificate for its status, or a
short reason why not.
"""
function check_certificate(lp::LP, sol::Solution)
    nv, nr = nvars(lp), nrows(lp)
    s = lp.sense === :max ? -1 : 1
    if sol.status === :optimal
        length(sol.x) == nv || return "x has the wrong length"
        length(sol.duals) == nr || return "duals have the wrong length"
        is_feasible(lp, sol.x) || return "x is not primal feasible"
        sol.objective == objective_value(lp, sol.x) || return "objective does not match x"
        y = s .* sol.duals                      # multipliers for the equivalent minimisation
        for (i, r) in enumerate(lp.rows)
            r.sense === :le && y[i] > 0 && return "dual of <= row $i has the wrong sign"
            r.sense === :ge && y[i] < 0 && return "dual of >= row $i has the wrong sign"
        end
        d = [s * lp.obj[j] for j in 1:nv] .- combine(lp, y)
        bound = s * lp.objconst + sum((y[i] * lp.rows[i].rhs for i in 1:nr); init = zero(Q))
        for j in 1:nv
            if d[j] > 0
                lp.lb[j] === nothing && return "positive reduced cost on variable $j without a lower bound"
                bound += d[j] * lp.lb[j]
            elseif d[j] < 0
                lp.ub[j] === nothing && return "negative reduced cost on variable $j without an upper bound"
                bound += d[j] * lp.ub[j]
            end
        end
        bound == s * sol.objective || return "dual bound $(fmt(bound)) differs from the primal objective"
        return nothing
    elseif sol.status === :infeasible
        conflicting_bounds(lp) != 0 && return nothing
        u = sol.farkas
        length(u) == nr || return "Farkas vector has the wrong length"
        for (i, r) in enumerate(lp.rows)
            r.sense === :le && u[i] < 0 && return "multiplier of <= row $i must be non-negative"
            r.sense === :ge && u[i] > 0 && return "multiplier of >= row $i must be non-positive"
        end
        w = combine(lp, u)
        lo = zero(Q)                            # smallest value of w'x over the bounds
        for j in 1:nv
            if w[j] > 0
                lp.lb[j] === nothing && return "combined row is unbounded below in variable $j"
                lo += w[j] * lp.lb[j]
            elseif w[j] < 0
                lp.ub[j] === nothing && return "combined row is unbounded below in variable $j"
                lo += w[j] * lp.ub[j]
            end
        end
        rhs = sum((u[i] * lp.rows[i].rhs for i in 1:nr); init = zero(Q))
        lo > rhs || return "combined row ($(fmt(lo)) <= $(fmt(rhs))) does not contradict the bounds"
        return nothing
    elseif sol.status === :unbounded
        is_feasible(lp, sol.x) || return "the starting point is not feasible"
        length(sol.ray) == nv || return "ray has the wrong length"
        r = sol.ray
        for j in 1:nv
            lp.lb[j] !== nothing && r[j] < 0 && return "ray decreases variable $j below its lower bound"
            lp.ub[j] !== nothing && r[j] > 0 && return "ray increases variable $j above its upper bound"
        end
        for (i, row) in enumerate(lp.rows)
            v = row_activity(row, r)
            row.sense === :le && v > 0 && return "ray violates <= row $i"
            row.sense === :ge && v < 0 && return "ray violates >= row $i"
            row.sense === :eq && v != 0 && return "ray violates = row $i"
        end
        gain = sum((s * lp.obj[j] * r[j] for j in 1:nv); init = zero(Q))
        gain < 0 || return "ray does not improve the objective"
        return nothing
    end
    "unknown status $(sol.status)"
end

"`true` when `sol`'s certificate holds for `lp` (see [`check_certificate`](@ref))."
verify(lp::LP, sol::Solution) = check_certificate(lp, sol) === nothing
