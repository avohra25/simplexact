# simplexact

Exact linear and mixed-integer programming in pure Julia. Every number is a
`Rational{BigInt}`, so statuses and optima are exact: no tolerances, no
"almost feasible", no wrong answers from rounding. Every answer also comes with
a certificate that an independent checker verifies against the original
problem.

Standard library only, no dependencies. Tested on Julia 1.11.

## What's in it

- **Two-phase primal simplex** on a dense exact tableau. Phase 1 uses
  artificial variables, drives the leftover ones out of the basis and keeps
  redundant equality rows harmless. Pivoting is Dantzig's rule with an
  automatic switch to Bland's rule while pivots stay degenerate, so it cannot
  cycle (Beale's classic cycling example is in the tests and examples).
- **Any bounds**: negative or fractional lower bounds, upper-only bounds, free
  variables, `<=`, `>=` and `=` rows, negative right-hand sides.
- **Certificates and an independent verifier** for all three outcomes:
  dual multipliers proving optimality, Farkas multipliers proving
  infeasibility, and a feasible point plus improving ray proving
  unboundedness.
- **Shadow prices and reduced costs** for optimal LPs.
- **Branch and bound** for mixed-integer programs with exact integrality,
  best-first or depth-first search, objective-lattice bound rounding, a node
  limit, and a proof tree that `verify_mip` re-checks leaf by leaf.
- **LP file format** (a CPLEX-style subset) with a parser that reports line
  numbers and a writer that round-trips exactly.
- **CLI** to solve, check and reformat LP files.

## Quick start

```julia
using Simplexact

lp = LP(:max)
x = add_variable!(lp, "x")                     # default bounds: 0 <= x < inf
y = add_variable!(lp, "y"; lb = -1, ub = "5/2")
set_objective!(lp, [x => 3, y => 2])
add_constraint!(lp, [x => 1, y => 1], :le, 4; name = "cap")
add_constraint!(lp, [x => 1, y => 3], :le, 6)
add_constraint!(lp, [x => 1], :le, 3)

sol = solve(lp)
sol.status          # :optimal
sol.objective       # 11//1
sol[x], sol[y]      # (3//1, 1//1)
sol.duals           # shadow price of each row (d objective / d rhs)
verify(lp, sol)     # true
```

Coefficients and bounds accept integers, rationals, strings such as `"7/3"` or
`"0.125"`, and floats (converted through their shortest decimal form, so `0.1`
becomes exactly `1//10`). Use `nothing` for an infinite bound.

## Certificates

`verify(lp, sol)` (or `check_certificate(lp, sol)`, which returns `nothing` or a
reason) uses only the original rows, bounds and objective, never the solver's
tableau, so a bug in the simplex code cannot make a wrong answer verify.

| status | certificate | what the verifier checks |
|---|---|---|
| `:optimal` | `x`, `duals` | `x` is feasible; the duals have the right signs and their dual bound (rows plus variable bounds) equals the objective, which is strong duality |
| `:infeasible` | `farkas` | the weighted sum of the rows is an inequality that no point within the variable bounds can satisfy |
| `:unbounded` | `x`, `ray` | `x` is feasible, moving along `ray` never violates a row or bound, and it strictly improves the objective |

```julia
lp = parse_lp("""
maximize
  x + y
subject to
  a: x + 2 y <= 4
  b: 3 x + y <= 6
  c: x + y >= 5
end
""")
sol = solve(lp)
sol.status          # :infeasible
sol.farkas          # [2//5, 1//5, -1]: 2/5 a + 1/5 b - c gives 0 <= -11/5
verify(lp, sol)     # true
```

## Mixed-integer programs

Mark variables with `integer = true` and call `solve_mip`:

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
that the two children of each branch cover every integer value, and re-checks
each leaf's LP certificate. A reported optimum is proven, not just found.

Other details:

- integer bounds are rounded inward before solving;
- when every variable in the objective is an integer variable, objective
  values lie on a lattice `objconst + g * Z`, and LP bounds are rounded up to
  that lattice before pruning, which closes nodes early without losing
  exactness;
- options: `strategy = :best_first | :depth_first`,
  `branching = :most_fractional | :first`, `node_limit`, and `rule` for the LP
  pivoting rule;
- a node limit returns `:node_limit` with the best incumbent so far and a
  proven bound;
- an unbounded relaxation is reported as `:relaxation_unbounded` instead of
  being guessed at.

## LP file format

```
\ comments start with a backslash
maximize
  profit: 3 x + 2 y + 1/2
subject to
  cap: x + y <= 4
  x + 3 y
     <= 6                  \ a row may span lines until its comparison
bounds
  x <= 3
  -1 <= y <= 5/2
  z free
general
  x
binary
  b
end
```

Variables default to `0 <= x < inf`. Numbers may be integers, decimals,
exponents or fractions like `2/3`; all are read exactly. Rows may have
variables on both sides. Ranged rows (`a <= expr <= b`) are rejected with a
clear message: write two rows. `read_lp(path)`, `parse_lp(text)` and
`write_lp(lp)` are the API; errors are `LPParseError`s with a line number.

## Command line

```
$ julia --project=. bin/simplexact.jl solve examples/diet.lp --decimal
status: optimal
objective: 25603/97  (≈ 263.948454)
simplex iterations: 15
  oats = 0
  milk = 491/97  (≈ 5.061856)
  bread = 12
  eggs = 0
  beans = 283/97  (≈ 2.917526)
certificate: verified

$ julia --project=. bin/simplexact.jl solve examples/knapsack.lp
status: optimal
objective: 171
nodes: 29, simplex iterations: 553
  ...
proof: verified (branch-and-bound tree checked exactly)
```

Commands: `solve`, `check` (summary) and `format` (canonical LP text). Solve
options: `--relax`, `--rule auto|bland|dantzig`,
`--strategy best_first|depth_first`, `--node-limit N`, `--decimal`, `--duals`
and `--no-verify`.

## Examples

`examples/` has a diet problem, a transportation problem, a 0/1 knapsack,
infeasible and unbounded LPs, Beale's cycling example, and `api.jl`, a small
production-planning MIP built from Julia:

```
julia --project=. examples/api.jl
```

## Tests

```
julia --project=. test/runtests.jl
```

About 5,300 checks, including:

- hundreds of random bounded LPs compared against brute-force vertex
  enumeration (every choice of active hyperplanes, solved exactly);
- hundreds more with open bounds, where every status (optimal, infeasible and
  unbounded all show up) must produce a certificate that verifies, and all
  three pivoting rules must agree;
- random MIPs compared against exhaustive enumeration of integer assignments,
  with every proof tree checked;
- tampered certificates and proof trees that the verifiers must reject;
- shadow prices checked against actually perturbing the right-hand side;
- LP-format round trips on random problems, parser error messages, and the CLI.

## Design notes

- Dense tableau, exact arithmetic: simple and fully trustworthy. It's meant for
  small and medium problems where exactness matters (teaching, verification,
  checking float solvers), not for large-scale speed.
- Standard-form conversion substitutes `x = l + x'`, `x = u - x'` or
  `x = x⁺ - x⁻` depending on the bounds and adds rows for finite two-sided
  boxes; certificates are mapped back to the original rows and variables, so
  users never see the internal form.
- Phase-1 duals are exactly a Farkas certificate, and the tableau's artificial
  columns hold the basis inverse, so all certificates come out of the solve for
  free.

## License

MIT, see [LICENSE](LICENSE).
