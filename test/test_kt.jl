K = 0.155 # [1/s]
T = 80.5 # [s]

start_time_second = 0.0
time_second_interval = 0.01
end_time_second = 200.0
time_list = start_time_second:time_second_interval:end_time_second

@testset "kt.jl" begin
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    kt_results = kt_simulate(K, T, time_list, δ_list)
end

@testset "kt_zigzag_test" begin
    target_δ_rad = 10.0 * π / 180.0
    target_ψ_rad_deviation = 10.0 * π / 180.0
    u_list, v_list, r_list, x_list, y_list, ψ_list, δ_list =
        kt_zigzag_test(K, T, time_list, target_δ_rad, target_ψ_rad_deviation)
end

@testset "kt least square method func" begin
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    u, v, r, x, y, Ψ, δ = kt_simulate(K, T, time_list, δ_list)

    data = ShipData(time_list, u, v, r, x, y, Ψ, δ, 0)
    K_est, T_est = estimate_kt_lsm(data)
    @test abs(K - K_est) < 1.0
    @test abs(T - T_est) < 10.0
end

@testset "kt bootstrap least square method func" begin
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    u, v, r, x, y, Ψ, δ = kt_simulate(K, T, time_list, δ_list)
    noize_dist = Normal(0.0, 0.0005)
    r_obs = r + rand(noize_dist, size(r))
    data = ShipData(time_list, u, v, r_obs, x, y, Ψ, δ, 0)
    one_sample_size = 5000
    K_est_samples, T_est_samples =
        estimate_kt_lsm_time_window_sampling(data, one_sample_size)
end

@testset "kt mcmc sampling func" begin
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    u, v, r, x, y, Ψ, δ = kt_simulate(K, T, time_list, δ_list)
    sampling_rate = 10
    time_obs = time_list[1:sampling_rate:end]
    noize_dist = Normal(0.0, 0.005)
    r_w_noize = r + rand(noize_dist, size(r))
    r_obs = r_w_noize[1:sampling_rate:end]
    δ_obs = δ[1:sampling_rate:end]
    data = ShipData(time_obs, 0, 0, r_obs, 0, 0, 0, δ_obs, 0)
    n_samples = 10
    n_chains = 1
    model = create_model_for_mcmc_sample_kt(data)
    chain = nuts_sampling_single_thread(model, n_samples, n_chains)
end
