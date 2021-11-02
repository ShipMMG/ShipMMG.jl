@testset "calc_position func" begin
    time_vec = [0.0, 2.0]
    u_vec = [1.0, 1.0]
    v_vec = [0.0, 0.0]
    r_vec = [0.0, 0.0]
    x, y, ψ = calc_position(time_vec, u_vec, v_vec, r_vec, x0 = 1.0, y0 = 1.0, ψ0 = 0.0)
    @test x[2] ≈ 3.0
    @test y[2] ≈ 1.0
    @test ψ[2] ≈ 0.0
end

@testset "draw_gif_result func" begin
    @testset "using KT" begin
        K_log = 0.155  # [1/s]
        T_log = 80.5  # [s]
        u0 = 10 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
        duration = 50  # [s]
        sampling = 1000
        time_list = range(0.0, stop = duration, length = sampling)
        Ts = 50.0
        δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
        kt_results = kt_simulate(K_log, T_log, time_list, δ_list)
        r, δ = kt_results
        u = u0 * ones(Float64, length(time_list))
        v = zeros(Float64, length(time_list))
        x, y, ψ = calc_position(time_list, u, v, r)

        test_result_file_name = "test_kt.gif"
        shape = [20, 5]
        draw_gif_result(time_list, x, y, ψ, shape, test_result_file_name)
        rm(test_result_file_name)
    end

end