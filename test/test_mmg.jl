basic_params, maneuvering_params = get_KVLCC2_L7_params()

@testset "mmg.jl KVLCC2_L7 turning" begin

    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]

    sampling = duration * 10
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        n_p_list,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    u, v, r, x, y, ψ, δ, n_p = mmg_results
    x_est, y_est, ψ_est = calc_position(time_list, u, v, r)
end

@testset "mmg_zigzag_test" begin
    target_δ_rad = 20.0 * π / 180.0
    target_ψ_rad_deviation = 20.0 * π / 180.0
    start_time_second = 0.00
    time_second_interval = 0.01
    end_time_second = 80.00
    time_list = start_time_second:time_second_interval:end_time_second
    n_const = 17.95  # [rps]
    n_p_list = n_const * ones(Float64, length(time_list))
    u_list, v_list, r_list, x_list, y_list, ψ_list, δ_list = mmg_3dof_zigzag_test(
        basic_params,
        maneuvering_params,
        time_list,
        n_p_list,
        target_δ_rad,
        target_ψ_rad_deviation,
    )
end

@testset "mmg mcmc sampling func" begin
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]

    sampling = duration * 10 + 1
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        n_p_list,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    u, v, r, δ, n_p = mmg_results
    sampling_rate = 10
    time_obs = time_list[1:sampling_rate:end]
    noize_dist = Normal(0.0, 0.01)
    u_obs = u + rand(noize_dist, size(u))
    v_obs = v + rand(noize_dist, size(v))
    r_obs = r + rand(noize_dist, size(r))
    x, y, ψ = calc_position(time_obs, u_obs, v_obs, r_obs)
    data = ShipData(
        time_obs,
        u_obs,
        v_obs,
        r_obs,
        x,
        y,
        ψ,
        δ[1:sampling_rate:end],
        n_p[1:sampling_rate:end],
    )
    n_samples = 10
    n_chains = 1
    model = create_model_for_mcmc_sample_mmg(
        data,
        basic_params,
        maneuvering_params.k_0,
        maneuvering_params.k_1,
        maneuvering_params.k_2,
    )
    chain = nuts_sampling_single_thread(model, n_samples, n_chains)
end
