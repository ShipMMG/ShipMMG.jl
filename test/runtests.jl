using ShipMMG
using Test

@testset "ShipMMG.jl" begin

    include("test_kt.jl")

    @testset "rotate_pos func" begin
        @test rotate_pos([4,0], deg2rad(90)) ≈ [0,4]
        @test rotate_pos([4,0], deg2rad(-90)) ≈ [0,-4]
    end

    @testset "square func" begin
        @test square(1, 0, [4,2], deg2rad(90)) ≈ [[0,0,2,2],[2,-2,-2,2]]
    end

    # @testset "draw_gif_result func" begin
    #     time, x, y, ψ, u, r, δ = kt_results
    # end

end
