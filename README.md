# simplexact

An exact linear-programming solver in pure Julia. Every number is a
`Rational{BigInt}`, so statuses and optima are exact: there are no tolerances,
no "almost feasible", and no wrong answers from rounding.

Standard library only, no dependencies.

## Status

Work in progress; this first commit has the model API and the two-phase primal
simplex method. More follows over the next few commits.

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

## Running the tests

```
julia --project test/runtests.jl
```

## License

MIT
