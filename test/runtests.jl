using Test
using Simplexact
include("helpers.jl")

@testset "Simplexact" begin
    include("test_model.jl")
    include("test_solve.jl")
    include("test_certificates.jl")
end
