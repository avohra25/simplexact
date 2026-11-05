# simplexact

An exact linear-programming solver in pure Julia. Every number is a
`Rational{BigInt}`, so statuses and optima are exact: there are no tolerances,
no "almost feasible", and no wrong answers from rounding.

Standard library only, no dependencies.

## Status

Work in progress. So far: the model API, the two-phase primal simplex method,
exact certificates with an independent verifier, and branch and bound for
mixed-integer programs with a checkable proof tree. A file format and CLI come
next.

```julia
using Simplexact

lp = LP(:max)
x = add_variable!(lp, "x")
y = add_variable!(lp, "y")
set_objective!(lp, [x => 3, y => 2])
add_constraint!(lp, [x => 1, y => 1], :le, 4)
add_constraint!(lp, [x => 1, y => 3], :le, 6)
add_constraint!(lp, [x => 1], :le, 3)

sol = solve(lp)
sol.status      # :optimal
sol.objective   # 11//1
sol[x], sol[y]  # (3//1, 1//1)
```

Variables take exact bounds (`lb`/`ub`, `nothing` for infinite, negative
lower bounds and free variables are fine). Rows can be `:le`, `:ge` or `:eq`.

## Certificates

Every answer comes with proof that can be checked without trusting the simplex
code. `verify(lp, sol)` re-checks it against the original problem using only
exact arithmetic:

- optimal: a feasible point plus dual multipliers (`sol.duals`, shadow prices,
  and `sol.reduced_costs`) whose dual bound equals the objective (strong duality);
- infeasible: Farkas multipliers (`sol.farkas`) whose combination of the rows
  contradicts the variable bounds;
- unbounded: a feasible point and a ray (`sol.ray`) that stays feasible while
  the objective improves.

```julia
sol = solve(lp)
verify(lp, sol)                 # true
check_certificate(lp, sol)      # nothing, or a string saying what is wrong
sol.duals                       # d objective / d rhs for each row
```

The test suite compares against brute-force vertex enumeration on hundreds of
random LPs and verifies certificates on thousands more, including unbounded
and infeasible ones.

## Mixed-integer programs

Mark variables with `integer = true` and call `solve_mip`. Integrality is
tested exactly, so there is no "0.9999999 is close enough to 1".

```julia
lp = LP(:max)
xs = [add_variable!(lp, "x$i"; ub = 1, integer = true) for i in 1:3]
set_objective!(lp, [xs[1] => 60, xs[2] => 100, xs[3] => 120])
add_constraint!(lp, [xs[1] => 10, xs[2] => 20, xs[3] => 30], :le, 50)

sol = solve_mip(lp)
sol.status, sol.objective   # (:optimal, 220//1)
verify_mip(lp, sol)         # true
```

`solve_mip` keeps the whole branch-and-bound tree as a proof. Each leaf holds
the LP solution for its box: a Farkas certificate (infeasible), an integral
optimum, or a dual bound that cannot beat the incumbent. `verify_mip`
recomputes every box from the root bounds and the branching decisions, checks
that the children of each branch cover every integer value, and re-checks each
leaf's LP certificate, so a reported optimum is proven, not just found.

Other details:

- integer bounds are rounded inward before solving;
- when every variable in the objective is integer, objective values live on a
  lattice (`objconst + g * Z`), and LP bounds are rounded up to it before
  pruning, which closes many nodes early without losing exactness;
- `strategy = :best_first | :depth_first`, `branching = :most_fractional | :first`;
- `node_limit` stops early with status `:node_limit`, the best incumbent so far
  and a proven bound;
- an unbounded relaxation is reported as `:relaxation_unbounded` rather than
  guessed at.

## Running the tests

```
julia --project test/runtests.jl
```

## License

MIT
