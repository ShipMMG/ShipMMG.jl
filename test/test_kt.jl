@testset "kt.jl" begin
    K_log = 0.155  # [1/s]
    T_log = 80.5  # [s]
    u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
    duration = 500  # [s]
    sampling = 10000
    time_list = range(0.0, stop = duration, length = sampling)
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    kt_results = kt_simulate(time_list, δ_list, K_log, T_log, u0)
end

@testset "kt_zigzag_test" begin
    K = 0.155
    T = 80.5
    target_δ_rad = 30.0 * π / 180.0
    target_ψ_rad_deviation = 10.0 * π / 180.0
    time_second_interval = 0.01
    end_time_second = 500.0
    u0 = 20 * (1852.0 / 3600) # [knot]
    δ_list, r_list, ψ_list = kt_zigzag_test(
        K,
        T,
        target_δ_rad,
        target_ψ_rad_deviation,
        time_second_interval,
        end_time_second,
        u0,
    )
end
