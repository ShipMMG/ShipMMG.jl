using ShipMMG
using Test

@testset "ShipMMG.jl" begin

    include("test_kt.jl")
    include("test_mmg.jl")
    include("test_draw.jl")

end
