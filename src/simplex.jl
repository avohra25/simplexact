# Two-phase primal simplex on a dense exact tableau.

mutable struct Tableau
    T::Matrix{Q}          # m x (n + m + 1): structural/slack columns, artificials, rhs
    basis::Vector{Int}
    n::Int
    m::Int
    iterations::Int
end

"Raw outcome of the standard-form solve."
struct StdResult
    status::Symbol                    # :optimal, :infeasible or :unbounded
    x::Vector{Q}                      # standard-form point (feasible when optimal or unbounded)
    y::Vector{Q}                      # row multipliers: optimal duals, or the Farkas vector when infeasible
    d::Vector{Q}                      # unbounded: improving ray in standard-form space
    basis::Vector{Int}
    iterations::Int
end

function pivot!(tab::Tableau, rc::Vector{Q}, r::Int, e::Int)
    T = tab.T
    ncols = size(T, 2)
    piv = T[r, e]
    if !isone(piv)
        for j in 1:ncols
            T[r, j] /= piv
        end
    end
    for i in 1:tab.m
        i == r && continue
        f = T[i, e]
        iszero(f) && continue
        for j in 1:ncols
            iszero(T[r, j]) || (T[i, j] -= f * T[r, j])
        end
    end
    f = rc[e]
    if !iszero(f)
        for j in 1:(ncols-1)
            iszero(T[r, j]) || (rc[j] -= f * T[r, j])
        end
    end
    tab.basis[r] = e
    tab.iterations += 1
    nothing
end

"Reduced costs `cost - cost_B' * T` for every non-rhs column."
function reduced_costs(tab::Tableau, cost::Vector{Q})
    rc = copy(cost)
    for i in 1:tab.m
        cb = cost[tab.basis[i]]
        iszero(cb) && continue
        for j in eachindex(rc)
            rc[j] -= cb * tab.T[i, j]
        end
    end
    rc
end

"""
Run primal simplex iterations on columns `1:ncand` until optimal or unbounded.
`rule` is `:dantzig`, `:bland` or `:auto` (Dantzig, switching to Bland's rule
while pivots keep being degenerate, which rules out cycling).
"""
function iterate!(tab::Tableau, rc::Vector{Q}, ncand::Int, rule::Symbol)
    T, m = tab.T, tab.m
    bland = rule === :bland
    stall = 0
    while true
        e = 0
        if bland
            for j in 1:ncand
                if rc[j] < 0
                    e = j
                    break
                end
            end
        else
            best = zero(Q)
            for j in 1:ncand
                if rc[j] < best
                    best = rc[j]
                    e = j
                end
            end
        end
        e == 0 && return (:optimal, 0)
        r = 0
        best = zero(Q)
        for i in 1:m
            a = T[i, e]
            a > 0 || continue
            ratio = T[i, end] / a
            if r == 0 || ratio < best || (ratio == best && tab.basis[i] < tab.basis[r])
                r, best = i, ratio
            end
        end
        r == 0 && return (:unbounded, e)
        pivot!(tab, rc, r, e)
        if rule === :auto
            if iszero(best)
                stall += 1
                stall > tab.m + tab.n && (bland = true)
            else
                stall = 0
                bland = false
            end
        end
    end
end

"Solve a standard-form problem with the two-phase primal simplex method."
function simplex_std(sf::StdForm; rule::Symbol = :auto)
    rule in (:auto, :dantzig, :bland) || throw(ArgumentError("unknown pivot rule :$rule"))
    m, n = size(sf.A)
    T = zeros(Q, m, n + m + 1)
    T[:, 1:n] = sf.A
    for i in 1:m
        T[i, n+i] = 1
    end
    T[:, end] = sf.b
    tab = Tableau(T, collect(n+1:n+m), n, m, 0)

    # Phase 1: minimise the sum of artificial variables.
    cost1 = zeros(Q, n + m)
    cost1[n+1:end] .= 1
    rc = reduced_costs(tab, cost1)
    iterate!(tab, rc, n + m, rule)
    infeas = sum((T[i, end] for i in 1:m if tab.basis[i] > n); init = zero(Q))
    if infeas > 0
        y = [1 - rc[n+i] for i in 1:m]
        return StdResult(:infeasible, Q[], y, Q[], copy(tab.basis), tab.iterations)
    end
    # Pivot remaining (zero-level) artificials out of the basis; a row with no
    # usable column is redundant and its artificial can stay at zero forever.
    for i in 1:m
        tab.basis[i] > n || continue
        for j in 1:n
            if !iszero(T[i, j])
                pivot!(tab, rc, i, j)
                break
            end
        end
    end

    # Phase 2: the real objective; artificials may not re-enter.
    cost2 = vcat(sf.c, zeros(Q, m))
    rc = reduced_costs(tab, cost2)
    status, e = iterate!(tab, rc, n, rule)
    x = zeros(Q, n)
    for i in 1:m
        tab.basis[i] <= n && (x[tab.basis[i]] = T[i, end])
    end
    if status === :unbounded
        d = zeros(Q, n)
        d[e] = 1
        for i in 1:m
            tab.basis[i] <= n && (d[tab.basis[i]] = -T[i, e])
        end
        return StdResult(:unbounded, x, Q[], d, copy(tab.basis), tab.iterations)
    end
    y = [-rc[n+i] for i in 1:m]
    StdResult(:optimal, x, y, Q[], copy(tab.basis), tab.iterations)
end
