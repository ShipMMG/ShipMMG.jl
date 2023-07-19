using ShipMMG
using Test
using Distributions

@testset "ShipMMG.jl" begin
    include("test_data.jl")
    include("test_kt.jl")
    include("test_mmg.jl")
end
