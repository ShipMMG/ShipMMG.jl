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

@testset "kt least square method func" begin
    K = 0.155  # [1/s]
    T = 80.5  # [s]
    duration = 500  # [s]
    sampling = 10000
    time_list = range(0.0, stop = duration, length = sampling)
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    r, δ = kt_simulate(K, T, time_list, δ_list)

    data = ShipData(time_list, 0, 0, r, 0, 0, 0, δ, 0)
    K_est, T_est = estimate_kt_lsm(data)
    @test abs(K - K_est) < 1.0
    @test abs(T - T_est) < 10.0
end
