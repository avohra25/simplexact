using Test
using Simplexact
include("helpers.jl")

@testset "Simplexact" begin
    include("test_model.jl")
    include("test_solve.jl")
    include("test_certificates.jl")
    include("test_mip.jl")
    include("test_lpformat.jl")
    include("test_mps.jl")
end
