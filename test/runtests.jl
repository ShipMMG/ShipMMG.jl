using JuliaShipMMG
using Test

@testset "JuliaShipMMG.jl" begin

    @testset "kt.jl" begin
        K_log = 0.155  # [1/s]
        T_log = 80.5  # [s]
        u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
        duration = 500  # [s]
        sampling = 10000
        time_list = range(0.0,stop=duration,length=sampling)
        Ts = 50.0
        δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
        @time sol=kt_simulate(time_list, δ_list, K_log, T_log, u0)
    end

end
