# simplexact

An exact linear-programming solver in pure Julia. Every number is a
`Rational{BigInt}`, so statuses and optima are exact: there are no tolerances,
no "almost feasible", and no wrong answers from rounding.

Standard library only, no dependencies.

## Status

Work in progress. So far: the model API, the two-phase primal simplex method,
and exact certificates with an independent verifier. Mixed-integer solving and
a file format/CLI come next.

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

## Running the tests

```
julia --project test/runtests.jl
```

## License

MIT
