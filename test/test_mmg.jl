basic_params, maneuvering_params = get_KVLCC2_L7_params()

@testset "mmg.jl KVLCC2_L7 turning" begin

    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]

    sampling = duration * 10
    time_list = range(0.00, stop=duration, length=sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        n_p_list,
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
    time_list = range(0.00, stop=duration, length=sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        n_p_list,
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

@testset "mmg mcmc multivariate sampling func" begin
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rps]

    sampling = duration * 10 + 1
    time_list = range(0.00, stop=duration, length=sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    n_p_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        n_p_list,
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

    parameter_width = 0.01
    μ = [
        0.022,    
        -0.040,   
        0.002,    
        0.011,    
        0.771,    
        -0.315,   
        0.083,    
        -1.607,   
        0.379,    
        -0.391,   
        0.008,    
        -0.137,   
        -0.049,   
        -0.030,   
        -0.294,   
        0.055,    
        -0.013    
    ]
    variance = parameter_width^2 / 3
    Σ = Diagonal(fill(variance, length(μ)))
    multivariate_prior_dist = MvNormal(μ, Σ)

    n_samples = 10
    n_chains = 1
    model = create_model_for_mcmc_sample_mmg(
        data,
        basic_params,
        maneuvering_params.k_0,
        maneuvering_params.k_1,
        maneuvering_params.k_2,
        multivariate_prior_dist,
        ρ=1025.0,
        σ_u_prior_dist=Uniform(0.00, 0.20),
        σ_v_prior_dist=Uniform(0.00, 0.20),
        σ_r_prior_dist=Uniform(0.00, 0.20),
    )
    chain = nuts_sampling_single_thread(model, n_samples, n_chains)
end

@testset "mpc for maneuvering variables" begin
    sampling_rate = 1
    duration = 100
    sampling = Int(duration * sampling_rate) + 1
    time_list = range(0.00, stop=duration, length=sampling)

    δ_list = 35 ./ 180.0 .* pi .* ones(Float64, sampling)
    n_p_list = 17.95 .* ones(Float64, sampling)

    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_list,
        n_p_list,
        u0=1.17,
        v0=0.0,
        r0=0.0,
    )
    u, v, r, x, y, Ψ, δ, n_p = mmg_results
    x1 = x .+ 3.0 .* cos.(Ψ)
    y1 = y .+ 3.0 .* sin.(Ψ)
    x2 = x .+ 3.0 .* cos.(Ψ .+ pi)
    y2 = y .+ 3.0 .* sin.(Ψ .+ pi)

    data = BowSternCoords(
        time = time_list,
        x1 = x1,
        y1 = y1,
        x2 = x2,
        y2 = y2,
    )

    t, u, v, r, x, y, Ψ, x1, y1, x2, y2, δ, n_p = mpc_for_maneuvering_variables(
        basic_params,
        maneuvering_params,
        data,
        u0=1.17,
        δ0=35/180.0*pi,
        n_p0=17.95,
        T_all = 50.0,
        T_step = 1.0,
        Np = 4,
    )    
end

@testset "mpc for external force" begin
    wind_force_and_moment_params = get_example_ship_wind_force_moment_params()
    sampling_rate = 4
    duration = 100
    sampling = Int(duration * sampling_rate) + 1
    time_list = range(0.00, stop=duration, length=sampling)

    δ_list = 35 ./ 180.0 .* pi .* ones(Float64, sampling)
    n_p_list = 17.95 .* ones(Float64, sampling)
    μ_U_W = 5.0
    μ_ψ_W = 1.0 * pi
    σ_U_W = 0.0
    σ_ψ_W = 0.0
    U_W_list = rand(Normal(μ_U_W, σ_U_W), length(time_list))
    ψ_W_list = rand(Normal(μ_ψ_W, σ_ψ_W), length(time_list))

    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        wind_force_and_moment_params,
        time_list,
        δ_list,
        n_p_list,
        U_W_list,
        ψ_W_list,
        u0=1.17,
        v0=0.0,
        r0=0.0,
    )
    u, v, r, x, y, Ψ, δ, n_p, X_wind, Y_wind, N_wind = mmg_results
    x1 = x .+ 3.0 .* cos.(Ψ)
    y1 = y .+ 3.0 .* sin.(Ψ)
    x2 = x .+ 3.0 .* cos.(Ψ .+ pi)
    y2 = y .+ 3.0 .* sin.(Ψ .+ pi)

    data = BowSternCoords(
        time = time_list,
        x1 = x1,
        y1 = y1,
        x2 = x2,
        y2 = y2,
    )

    t, u, v, r, x, y, Ψ, x1, y1, x2, y2, δ, n_p, X_F, Y_F, N_F= mpc_for_external_force(
        basic_params,
        maneuvering_params,
        data,
        δ,
        n_p,
        u0=1.17,
        T_all = 50.0,
        T_step = 0.25,
        Np = 8,
    )
end