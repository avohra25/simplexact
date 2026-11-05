"""
Result of [`solve`](@ref).

`status` is `:optimal`, `:infeasible` or `:unbounded`. For an optimal solve,
`objective` and `x` hold the exact optimum; for an unbounded one `ray` is a
direction along which the objective improves forever.
"""
struct Solution
    status::Symbol
    objective::Union{Nothing,Q}
    x::Vector{Q}
    ray::Vector{Q}
    iterations::Int
end

Base.getindex(s::Solution, v::Var) = s.x[v.index]

"""
    solve(lp; rule = :auto) -> Solution

Solve the continuous relaxation of `lp` exactly with a two-phase primal simplex
method over `Rational{BigInt}`.
"""
function solve(lp::LP; rule::Symbol = :auto)
    conflicting_bounds(lp) == 0 || return Solution(:infeasible, nothing, Q[], Q[], 0)
    sf = to_standard(lp)
    res = simplex_std(sf; rule)
    if res.status === :infeasible
        Solution(:infeasible, nothing, Q[], Q[], res.iterations)
    elseif res.status === :unbounded
        Solution(:unbounded, nothing, Q[], original_direction(sf, res.d), res.iterations)
    else
        x = original_point(sf, res.x)
        Solution(:optimal, objective_value(lp, x), x, Q[], res.iterations)
    end
end
