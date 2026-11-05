using Test
using Simplexact

@testset "Simplexact" begin
    include("test_model.jl")
    include("test_solve.jl")
end
