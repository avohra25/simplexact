# Branch and bound for mixed-integer programs, with a checkable optimality proof.

"""
One node of a branch-and-bound proof tree. The box of a node is never stored:
it follows from the root bounds and the branching decisions above it, so a
verifier can recompute it instead of trusting the solver.

`kind` is `:branch` (children `down`: `x[var] <= split`, `up`: `x[var] >= split + 1`),
`:infeasible`, `:cutoff` (LP bound cannot beat the incumbent) or `:integral`.
Leaves carry the LP [`Solution`](@ref) for their box.
"""
mutable struct ProofNode
    kind::Symbol
    solution::Union{Nothing,Solution}
    var::Int
    split::Q
    down::Union{Nothing,ProofNode}
    up::Union{Nothing,ProofNode}
end
ProofNode(kind::Symbol, sol = nothing) = ProofNode(kind, sol, 0, zero(Q), nothing, nothing)

"""
Result of [`solve_mip`](@ref). `status` is `:optimal`, `:infeasible`,
`:node_limit` (best known answer so far in `x`/`objective`, proven bound in
`bound`) or `:relaxation_unbounded`. `proof` is a [`ProofNode`](@ref) tree that
[`verify_mip`](@ref) can check.
"""
struct MIPSolution
    status::Symbol
    objective::Union{Nothing,Q}
    x::Vector{Q}
    bound::Union{Nothing,Q}
    nodes::Int
    lp_iterations::Int
    proof::Union{Nothing,ProofNode}
end

Base.getindex(s::MIPSolution, v::Var) = s.x[v.index]

isint(q::Rational) = isone(denominator(q))

"A copy of `lp` with different variable bounds (rows and objective are shared)."
with_bounds(lp::LP, lb::Vector{Bound}, ub::Vector{Bound}) =
    LP(lp.sense, lp.names, lb, ub, lp.integer, lp.obj, lp.objconst, lp.rows)

"Root bounds: integer variables get their bounds rounded inward."
function tightened_bounds(lp::LP)
    lb, ub = copy(lp.lb), copy(lp.ub)
    for j in 1:nvars(lp)
        lp.integer[j] || continue
        lb[j] === nothing || (lb[j] = Q(ceil(BigInt, lb[j])))
        ub[j] === nothing || (ub[j] = Q(floor(BigInt, ub[j])))
    end
    lb, ub
end

"""
Spacing of attainable objective values: if every variable with a non-zero
objective coefficient is an integer variable, objective values all lie in
`objconst + step * Z` and LP bounds can be rounded up to that lattice.
Returns `nothing` when there is no such structure.
"""
function objective_step(lp::LP)
    nz = [j for j in 1:nvars(lp) if !iszero(lp.obj[j])]
    (isempty(nz) || !all(j -> lp.integer[j], nz)) && return nothing
    den = reduce(lcm, (denominator(lp.obj[j]) for j in nz))
    g = reduce(gcd, (numerator(lp.obj[j] * den) for j in nz))
    Q(g, den)
end

"Smallest attainable value of the (min-sense) objective that is at least `z`."
function lattice_ceil(lp::LP, step::Union{Nothing,Q}, z::Q)
    step === nothing && return z
    o = (lp.sense === :max ? -1 : 1) * lp.objconst
    o + step * Q(ceil(BigInt, (z - o) / step))
end

"""
    solve_mip(lp; rule = :auto, branching = :most_fractional, strategy = :best_first,
              node_limit = 100_000) -> MIPSolution

Exact branch and bound over the LP relaxation. Because the arithmetic is
rational, integrality is tested exactly (no tolerance). `branching` is
`:most_fractional` or `:first`; `strategy` is `:best_first` (lowest bound
next) or `:depth_first`.
"""
function solve_mip(lp::LP; rule::Symbol = :auto, branching::Symbol = :most_fractional,
    strategy::Symbol = :best_first, node_limit::Int = 100_000)
    branching in (:most_fractional, :first) || throw(ArgumentError("unknown branching rule :$branching"))
    strategy in (:best_first, :depth_first) || throw(ArgumentError("unknown strategy :$strategy"))
    s = lp.sense === :max ? -1 : 1
    step = objective_step(lp)
    ints = [j for j in 1:nvars(lp) if lp.integer[j]]
    nodes = 0
    iters = 0
    incumbent = nothing                    # (min-sense value, x)
    unbounded = false

    function evaluate(lb, ub)
        sol = solve(with_bounds(lp, lb, ub); rule)
        nodes += 1
        iters += sol.iterations
        sol
    end
    cutoff(sol) = incumbent !== nothing && lattice_ceil(lp, step, s * sol.objective) >= incumbent[1]
    function fractional(sol)
        cand = [j for j in ints if !isint(sol.x[j])]
        isempty(cand) && return 0
        branching === :first && return first(cand)
        argmin(j -> abs(sol.x[j] - floor(sol.x[j]) - 1 // 2), cand)
    end

    # Open nodes: (id, depth, lb, ub, solution, proof node).
    open = Tuple{Int,Q,Vector{Bound},Vector{Bound},Solution,ProofNode}[]
    nextid = 0
    function admit(lb, ub, sol)::ProofNode
        if sol.status === :infeasible
            return ProofNode(:infeasible, sol)
        elseif sol.status === :unbounded
            unbounded = true
            return ProofNode(:infeasible, sol)
        end
        if fractional(sol) == 0
            v = s * sol.objective
            if incumbent === nothing || v < incumbent[1]
                incumbent = (v, sol.x)
            end
            return ProofNode(:integral, sol)
        end
        node = ProofNode(:cutoff, sol)         # becomes :branch if it is ever expanded
        nextid += 1
        push!(open, (nextid, s * sol.objective, lb, ub, sol, node))
        node
    end

    lb0, ub0 = tightened_bounds(lp)
    root = admit(lb0, ub0, evaluate(lb0, ub0))
    while !unbounded && !isempty(open) && nodes < node_limit
        k = if strategy === :depth_first
            argmax(i -> open[i][1], eachindex(open))
        else
            argmin(i -> (lattice_ceil(lp, step, open[i][2]), -open[i][1]), eachindex(open))
        end
        (_, _, lb, ub, sol, node) = open[k]
        deleteat!(open, k)
        cutoff(sol) && continue                # stays a :cutoff leaf
        j = fractional(sol)
        f = Q(floor(BigInt, sol.x[j]))
        dub = copy(ub); dub[j] = f
        ulb = copy(lb); ulb[j] = f + 1
        node.kind = :branch
        node.var, node.split = j, f
        node.down = admit(lb, dub, evaluate(lb, dub))
        node.up = admit(ulb, ub, evaluate(ulb, ub))
    end

    if unbounded
        return MIPSolution(:relaxation_unbounded, nothing, Q[], nothing, nodes, iters, nothing)
    end
    open_bound = isempty(open) ? nothing : minimum(lattice_ceil(lp, step, t[2]) for t in open)
    if !isempty(open)                          # stopped early: the open nodes are unfinished
        for t in open
            t[6].kind = :open
        end
        best = incumbent === nothing ? nothing : incumbent[1]
        bound = best === nothing ? open_bound : min(open_bound, best)
        return MIPSolution(:node_limit, best === nothing ? nothing : s * best,
            incumbent === nothing ? Q[] : incumbent[2], s * bound, nodes, iters, nothing)
    end
    if incumbent === nothing
        return MIPSolution(:infeasible, nothing, Q[], nothing, nodes, iters, root)
    end
    MIPSolution(:optimal, s * incumbent[1], incumbent[2], s * incumbent[1], nodes, iters, root)
end

"""
    check_mip_certificate(lp, sol) -> Union{Nothing,String}

Check the proof tree of an optimal or infeasible [`MIPSolution`](@ref) with exact
arithmetic: the incumbent must be feasible and integral, the tree must cover
every integer assignment, and every leaf must be justified by a valid LP
certificate (see [`check_certificate`](@ref)).
"""
function check_mip_certificate(lp::LP, sol::MIPSolution)
    sol.status in (:optimal, :infeasible) || return "status $(sol.status) carries no proof"
    sol.proof === nothing && return "missing proof tree"
    s = lp.sense === :max ? -1 : 1
    step = objective_step(lp)
    zinc = nothing
    if sol.status === :optimal
        is_feasible(lp, sol.x) || return "incumbent is not feasible"
        all(j -> !lp.integer[j] || isint(sol.x[j]), 1:nvars(lp)) || return "incumbent is not integral"
        sol.objective == objective_value(lp, sol.x) || return "objective does not match the incumbent"
        zinc = s * sol.objective
    end
    function walk(node::ProofNode, lb, ub)
        node.kind === :open && return "proof contains an unfinished node"
        if node.kind === :branch
            j = node.var
            (1 <= j <= nvars(lp) && lp.integer[j]) || return "branching on a non-integer variable"
            isint(node.split) || return "non-integral split value"
            (node.down === nothing || node.up === nothing) && return "branch without children"
            dub = copy(ub); dub[j] = ub[j] === nothing ? node.split : min(ub[j], node.split)
            ulb = copy(lb); ulb[j] = lb[j] === nothing ? node.split + 1 : max(lb[j], node.split + 1)
            r = walk(node.down, lb, dub)
            r === nothing || return r
            return walk(node.up, ulb, ub)
        end
        node.solution === nothing && return "leaf without an LP solution"
        sub = with_bounds(lp, lb, ub)
        reason = check_certificate(sub, node.solution)
        reason === nothing || return "leaf certificate invalid: $reason"
        if node.kind === :infeasible
            return node.solution.status === :infeasible ? nothing : "leaf marked infeasible but is not"
        end
        node.solution.status === :optimal || return "leaf is not an optimal LP"
        v = s * node.solution.objective
        if node.kind === :cutoff
            zinc === nothing && return "cutoff leaf in a proof without an incumbent"
            lattice_ceil(lp, step, v) >= zinc || return "cutoff leaf could still beat the incumbent"
        elseif node.kind === :integral
            all(j -> !lp.integer[j] || isint(node.solution.x[j]), 1:nvars(lp)) || return "integral leaf is fractional"
            zinc === nothing && return "integral leaf in a proof without an incumbent"
            v >= zinc || return "integral leaf beats the claimed optimum"
        else
            return "unknown node kind $(node.kind)"
        end
        nothing
    end
    walk(sol.proof, tightened_bounds(lp)...)
end

"`true` when the proof tree of `sol` checks out for `lp`."
verify_mip(lp::LP, sol::MIPSolution) = check_mip_certificate(lp, sol) === nothing
