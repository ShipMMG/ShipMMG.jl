basic_params, maneuvering_params = get_KVLCC2_L7_params()
structure_params = get_structure_params()

@testset "mmg.jl KVLCC2_L7 turning" begin

    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]
    μ_u_wind = 8.0 # [m/s]
    μ_ψ_wind = 0.0 # [deg]  
    σ_u_wind = 1.0 # [m/s]
    σ_ψ_wind = 0.1 # [deg]

    sampling = duration * 10
    time_list = range(0.00, stop=duration, length=sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    u_wind_list = rand(Normal(μ_u_wind,σ_u_wind),length(time_list))
    ψ_wind_list = rand(Normal(μ_ψ_wind,σ_ψ_wind),length(time_list))

    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        structure_params,
        time_list,
        δ_rad_list,
        n_p_list,
        u_wind_list,
        ψ_wind_list,
        u0=2.29 * 0.512,
        v0=0.0,
        r0=0.0,
    )
    u, v, r, x, y, ψ, δ, n_p = mmg_results
    x_est, y_est, ψ_est = calc_position(time_list, u, v, r)
end

@testset "mmg_zigzag_test" begin
    target_δ_rad = 20.0 * π / 180.0
    target_ψ_rad_deviation = 20.0 * π / 180.0
    μ_u_wind = 8.0 # [m/s]
    μ_ψ_wind = 0.0 # [deg]  
    σ_u_wind = 1.0 # [m/s]
    σ_ψ_wind = 0.1 # [deg]
    start_time_second = 0.00
    time_second_interval = 0.01
    end_time_second = 80.00
    time_list = start_time_second:time_second_interval:end_time_second
    n_const = 17.95  # [rps]
    n_p_list = n_const * ones(Float64, length(time_list))
    u_wind_list = rand(Normal(μ_u_wind,σ_u_wind),length(time_list))
    ψ_wind_list = rand(Normal(μ_ψ_wind,σ_ψ_wind),length(time_list))
    u_list, v_list, r_list, x_list, y_list, ψ_list, δ_list = mmg_3dof_zigzag_test(
        basic_params,
        maneuvering_params,
        structure_params,
        time_list,
        n_p_list,
        target_δ_rad,
        target_ψ_rad_deviation,
        u_wind_list,
        ψ_wind_list,
    )
end

@testset "mmg mcmc sampling func" begin
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]
    μ_u_wind = 8.0 # [m/s]
    μ_ψ_wind = 0.0 # [deg]   
    σ_u_wind = 1.0 # [m/s]
    σ_ψ_wind = 0.1 # [deg]

    sampling = duration * 10 + 1
    time_list = range(0.00, stop=duration, length=sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    u_wind_list = rand(Normal(μ_u_wind,σ_u_wind),length(time_list))
    ψ_wind_list = rand(Normal(μ_ψ_wind,σ_ψ_wind),length(time_list))
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        structure_params,
        time_list,
        δ_rad_list,
        n_p_list,
        u_wind_list,
        ψ_wind_list,
        u0=2.29 * 0.512,
        v0=0.0,
        r0=0.0,
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
        δ_rad_list[1:sampling_rate:end],
        n_p_list[1:sampling_rate:end],
    )

    parameter_width = 0.001
    n_samples = 10
    n_chains = 1
    model = create_model_for_mcmc_sample_mmg(
        data,
        basic_params,
        maneuvering_params.k_0,
        maneuvering_params.k_1,
        maneuvering_params.k_2,
        ρ=1025.0,
        σ_u_prior_dist=Uniform(0.00, 0.20),
        σ_v_prior_dist=Uniform(0.00, 0.20),
        σ_r_prior_dist=Uniform(0.00, 0.20),
        R_0_dash_prior_dist=Uniform(0.022 - parameter_width, 0.022 + parameter_width),
        X_vv_dash_prior_dist=Uniform(-0.04 - parameter_width, -0.04 + parameter_width),
        X_vr_dash_prior_dist=Uniform(0.002 - parameter_width, 0.002 + parameter_width),
        X_rr_dash_prior_dist=Uniform(0.011 - parameter_width, 0.011 + parameter_width),
        X_vvvv_dash_prior_dist=Uniform(0.771 - parameter_width, 0.771 + parameter_width),
        Y_v_dash_prior_dist=Uniform(-0.315 - parameter_width, -0.315 + parameter_width),
        Y_r_dash_prior_dist=Uniform(0.083 - parameter_width, 0.083 + parameter_width),
        Y_vvv_dash_prior_dist=Uniform(-1.607 - parameter_width, -1.607 + parameter_width),
        Y_vvr_dash_prior_dist=Uniform(0.379 - parameter_width, 0.379 + parameter_width),
        Y_vrr_dash_prior_dist=Uniform(-0.391 - parameter_width, -0.391 + parameter_width),
        Y_rrr_dash_prior_dist=Uniform(0.008 - parameter_width, 0.008 + parameter_width),
        N_v_dash_prior_dist=Uniform(-0.137 - parameter_width, -0.137 + parameter_width),
        N_r_dash_prior_dist=Uniform(-0.049 - parameter_width, -0.049 + parameter_width),
        N_vvv_dash_prior_dist=Uniform(-0.03 - parameter_width, -0.03 + parameter_width),
        N_vvr_dash_prior_dist=Uniform(-0.294 - parameter_width, -0.294 + parameter_width),
        N_vrr_dash_prior_dist=Uniform(0.055 - parameter_width, 0.055 + parameter_width),
        N_rrr_dash_prior_dist=Uniform(-0.013 - parameter_width, -0.013 + parameter_width),
    )
    chain = nuts_sampling_single_thread(model, n_samples, n_chains)
end
