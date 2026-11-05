# Test helpers: a seeded generator of small LPs and an exact brute-force oracle.

using Random

"Solve the square system `A x = b` exactly; `nothing` if singular."
function gauss(A::Matrix{Q}, b::Vector{Q})
    n = length(b)
    M = hcat(copy(A), b)
    for c in 1:n
        p = findfirst(i -> !iszero(M[i, c]), c:n)
        p === nothing && return nothing
        p += c - 1
        M[c, :], M[p, :] = M[p, :], M[c, :]
        M[c, :] ./= M[c, c]
        for i in 1:n
            i != c && !iszero(M[i, c]) && (M[i, :] .-= M[i, c] .* M[c, :])
        end
    end
    M[:, end]
end

function subsets(n::Int, k::Int)
    k == 0 && return [Int[]]
    res = Vector{Vector{Int}}()
    function rec(start, cur)
        length(cur) == k && (push!(res, copy(cur)); return)
        for i in start:n
            push!(cur, i); rec(i + 1, cur); pop!(cur)
        end
    end
    rec(1, Int[])
    res
end

"""
Brute-force optimum of an LP whose variables all have finite bounds: every
vertex of the (bounded) polytope lies on `n` linearly independent active
hyperplanes, so enumerate all choices. Returns `(status, objective)`.
"""
function brute_force(lp::LP)
    n = nvars(lp)
    planes = Tuple{Vector{Q},Q}[]
    for r in lp.rows
        push!(planes, ([get(r.coeffs, j, zero(Q)) for j in 1:n], r.rhs))
    end
    for j in 1:n
        e = zeros(Q, n); e[j] = 1
        push!(planes, (e, lp.lb[j]))
        push!(planes, (e, lp.ub[j]))
    end
    best = nothing
    for S in subsets(length(planes), n)
        A = Matrix{Q}(undef, n, n)
        for (k, i) in enumerate(S)
            A[k, :] = planes[i][1]
        end
        x = gauss(A, Q[planes[i][2] for i in S])
        x === nothing && continue
        is_feasible(lp, x) || continue
        v = objective_value(lp, x)
        if best === nothing || (lp.sense === :min ? v < best : v > best)
            best = v
        end
    end
    best === nothing ? (:infeasible, nothing) : (:optimal, best)
end

"A random small LP. With `boxed` every variable has finite bounds."
function random_lp(rng; boxed::Bool = true, maxn = 3, maxm = 4)
    n = rand(rng, 1:maxn)
    m = rand(rng, 0:maxm)
    lp = LP(rand(rng, (:min, :max)))
    vs = Var[]
    for j in 1:n
        lo = rand(rng, -3:2)
        hi = lo + rand(rng, 0:5)
        lb = boxed || rand(rng) < 0.6 ? lo : nothing
        ub = boxed || rand(rng) < 0.6 ? hi : nothing
        push!(vs, add_variable!(lp, "v$j"; lb, ub))
    end
    set_objective!(lp, [vs[j] => rand(rng, -4:4) for j in 1:n]; constant = rand(rng, -3:3))
    for _ in 1:m
        add_constraint!(lp, [vs[j] => rand(rng, -3:3) for j in 1:n],
            rand(rng, (:le, :ge, :eq, :le, :ge)), rand(rng, -6:8))
    end
    lp
end
