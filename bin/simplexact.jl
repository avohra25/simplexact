#!/usr/bin/env julia
# Usage: julia --project=. bin/simplexact.jl solve examples/knapsack.lp
pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Simplexact
exit(Simplexact.main(ARGS))
