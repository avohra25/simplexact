# Conversion of an LP with arbitrary bounds into standard form
#     min c'x   s.t.  A x = b,  x >= 0,  b >= 0.

"""
Standard form of an [`LP`](@ref) plus the data needed to map solutions, rays
and multipliers back to the original variables and rows.
"""
struct StdForm
    A::Matrix{Q}
    b::Vector{Q}
    c::Vector{Q}
    c0::Q                                  # constant term of the (min-sense) objective
    nstruct::Int                           # number of structural columns (before slacks)
    offset::Vector{Q}                      # x_j = offset[j] + sum(k * x'[col] for (col, k) in cols[j])
    cols::Vector{Vector{Tuple{Int,Q}}}
    flip::Vector{Int}                      # -1 where a row was negated to make b >= 0
    nrows_orig::Int                        # rows 1..nrows_orig are the LP's own rows, the rest are bound rows
    sign::Int                              # +1 for min, -1 for max
end

"A pair of contradictory bounds `lb > ub`, if any (such an LP is trivially infeasible)."
function conflicting_bounds(lp::LP)
    for j in 1:nvars(lp)
        l, u = lp.lb[j], lp.ub[j]
        l !== nothing && u !== nothing && l > u && return j
    end
    return 0
end

function to_standard(lp::LP)
    conflicting_bounds(lp) == 0 || throw(ArgumentError("contradictory bounds"))
    s = lp.sense === :max ? -1 : 1
    nv = nvars(lp)
    cols = [Tuple{Int,Q}[] for _ in 1:nv]
    offset = zeros(Q, nv)
    ncol = 0
    bound_rows = Tuple{Int,Q}[]            # (std column, upper limit of the shifted variable)
    for j in 1:nv
        l, u = lp.lb[j], lp.ub[j]
        if l !== nothing                    # x = l + x'
            ncol += 1
            push!(cols[j], (ncol, one(Q)))
            offset[j] = l
            u !== nothing && push!(bound_rows, (ncol, u - l))
        elseif u !== nothing                # x = u - x'
            ncol += 1
            push!(cols[j], (ncol, -one(Q)))
            offset[j] = u
        else                                # x = x+ - x-
            ncol += 2
            push!(cols[j], (ncol - 1, one(Q)))
            push!(cols[j], (ncol, -one(Q)))
        end
    end
    nr = nrows(lp)
    m = nr + length(bound_rows)
    nslack = count(r -> r.sense !== :eq, lp.rows) + length(bound_rows)
    n = ncol + nslack
    A = zeros(Q, m, n)
    b = zeros(Q, m)
    flip = ones(Int, m)
    sl = ncol
    for (i, r) in enumerate(lp.rows)
        rhs = r.rhs
        for (j, a) in r.coeffs
            rhs -= a * offset[j]
            for (col, k) in cols[j]
                A[i, col] += a * k
            end
        end
        if r.sense === :le
            sl += 1; A[i, sl] = 1
        elseif r.sense === :ge
            sl += 1; A[i, sl] = -1
        end
        if rhs < 0
            A[i, :] .*= -1
            rhs = -rhs
            flip[i] = -1
        end
        b[i] = rhs
    end
    for (k, (col, w)) in enumerate(bound_rows)
        i = nr + k
        A[i, col] = 1
        sl += 1
        A[i, sl] = 1
        b[i] = w
    end
    c = zeros(Q, n)
    c0 = s * lp.objconst
    for j in 1:nv
        cj = s * lp.obj[j]
        c0 += cj * offset[j]
        for (col, k) in cols[j]
            c[col] += cj * k
        end
    end
    StdForm(A, b, c, c0, ncol, offset, cols, flip, nr, s)
end

"Map a standard-form point back to the original variables."
function original_point(sf::StdForm, xs::AbstractVector)
    [sf.offset[j] + sum((k * xs[col] for (col, k) in sf.cols[j]); init = zero(Q)) for j in eachindex(sf.cols)]
end

"Map a standard-form direction back to the original variables (no offsets)."
function original_direction(sf::StdForm, ds::AbstractVector)
    [sum((k * ds[col] for (col, k) in sf.cols[j]); init = zero(Q)) for j in eachindex(sf.cols)]
end
