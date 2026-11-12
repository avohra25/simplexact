"""
Simplexact: an exact linear-programming solver.

All arithmetic is done in `Rational{BigInt}`, so statuses and optima are exact
and never depend on floating-point tolerances.
"""
module Simplexact

export LP, Var, Solution, Q, add_variable!, add_constraint!, set_objective!, solve, verify, check_certificate,
    MIPSolution, ProofNode, solve_mip, verify_mip, check_mip_certificate,
    parse_lp, read_lp, write_lp, LPParseError,
    parse_mps, read_mps, write_mps, MPSParseError, read_model,
    parse_q, tonum, fmt, decimal, objective_value, is_feasible, nvars, nrows

include("model.jl")
include("stdform.jl")
include("simplex.jl")
include("solve.jl")
include("verify.jl")
include("mip.jl")
include("lpformat.jl")
include("mpsformat.jl")
include("cli.jl")

end # module
