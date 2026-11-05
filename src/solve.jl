"""
Result of [`solve`](@ref).

`status` is `:optimal`, `:infeasible` or `:unbounded`. Besides the answer, a
solution carries a certificate that [`verify`](@ref) checks independently of
the simplex code, using only exact arithmetic on the original problem:

* `:optimal`: `x` and `objective`, plus `duals` (one shadow price per row,
  i.e. d objective / d rhs) and `reduced_costs` (one per variable). A feasible
  `x` and a matching dual bound prove optimality (strong duality).
* `:infeasible`: `farkas`, one multiplier per row. Combining the rows with these
  multipliers gives an inequality that no point inside the variable bounds can
  satisfy. When the bounds themselves contradict each other `farkas` is empty.
* `:unbounded`: a feasible point `x` and a `ray` along which the objective
  improves forever without leaving the feasible region.
"""
struct Solution
    status::Symbol
    objective::Union{Nothing,Q}
    x::Vector{Q}
    ray::Vector{Q}
    duals::Vector{Q}
    reduced_costs::Vector{Q}
    farkas::Vector{Q}
    iterations::Int
end

Solution(status::Symbol; objective = nothing, x = Q[], ray = Q[], duals = Q[],
    reduced_costs = Q[], farkas = Q[], iterations = 0) =
    Solution(status, objective, x, ray, duals, reduced_costs, farkas, iterations)

Base.getindex(s::Solution, v::Var) = s.x[v.index]

"""
    solve(lp; rule = :auto) -> Solution

Solve the continuous relaxation of `lp` exactly with a two-phase primal simplex
method over `Rational{BigInt}`. `rule` picks the pivoting rule: `:dantzig`,
`:bland`, or `:auto` (Dantzig with a Bland fallback on degenerate stalls).
"""
function solve(lp::LP; rule::Symbol = :auto)
    conflicting_bounds(lp) == 0 || return Solution(:infeasible)
    sf = to_standard(lp)
    res = simplex_std(sf; rule)
    nr = sf.nrows_orig
    if res.status === :infeasible
        # Standard-form Farkas vector y has A'y <= 0 and b'y > 0; flip signs back
        # to original rows and negate so that <= rows get non-negative weights.
        u = [-sf.flip[i] * res.y[i] for i in 1:nr]
        Solution(:infeasible; farkas = u, iterations = res.iterations)
    elseif res.status === :unbounded
        Solution(:unbounded; x = original_point(sf, res.x), ray = original_direction(sf, res.d),
            iterations = res.iterations)
    else
        x = original_point(sf, res.x)
        ymin = [sf.flip[i] * res.y[i] for i in 1:nr]
        duals = sf.sign .* ymin
        rcs = lp.obj .- combine(lp, duals)
        Solution(:optimal; objective = objective_value(lp, x), x, duals, reduced_costs = rcs,
            iterations = res.iterations)
    end
end
