# Build and solve a small production-planning MIP from Julia.
# Run with: julia --project=. examples/api.jl
using Simplexact

lp = LP(:max)
chairs = add_variable!(lp, "chairs"; integer = true)
tables = add_variable!(lp, "tables"; integer = true)
overtime = add_variable!(lp, "overtime"; ub = 10)         # hours, continuous

set_objective!(lp, [chairs => 45, tables => 80, overtime => -30])
add_constraint!(lp, [chairs => 5, tables => 20], :le, 410; name = "wood")
add_constraint!(lp, [chairs => 10, tables => 15, overtime => -1], :le, 457; name = "labour")

relax = solve(lp)
println("LP relaxation: ", relax.status, ", objective ", fmt(relax.objective), " (verified: ", verify(lp, relax), ")")
println("  shadow prices: ", join(["$(r.name) = $(fmt(d))" for (r, d) in zip(lp.rows, relax.duals)], ", "))

sol = solve_mip(lp)
println("MIP: ", sol.status, ", objective ", fmt(sol.objective), " after ", sol.nodes, " LP solves")
println("  chairs = ", fmt(sol[chairs]), ", tables = ", fmt(sol[tables]), ", overtime = ", fmt(sol[overtime]))
println("  proof verified: ", verify_mip(lp, sol))
