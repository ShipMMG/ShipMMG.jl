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

    @testset "draw_gif_result func" begin
        @testset "using KT" begin
            K_log = 0.155  # [1/s]
            T_log = 80.5  # [s]
            u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
            duration = 50  # [s]
            sampling = 1000
            time_list = range(0.0,stop=duration,length=sampling)
            Ts = 50.0
            δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
            kt_results=kt_simulate(time_list, δ_list, K_log, T_log, u0)
            time, x, y, ψ, u, r, δ = kt_results

            test_result_file_name = "test_kt.gif"
            draw_gif_result(time, x, y, ψ,test_result_file_name)
            rm(test_result_file_name)
        end
        
    end

end
