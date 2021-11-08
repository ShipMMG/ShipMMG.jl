@with_kw struct ShipData{Tt,Tu,Tv,Tr,Tx,Ty,Tψ,Tδ,Tnpm}
    time::Tt
    u::Tu
    v::Tv
    r::Tr
    x::Tx
    y::Ty
    ψ::Tψ
    δ::Tδ
    npm::Tnpm
end

function calc_position(time_vec, u_vec, v_vec, r_vec; x0 = 0.0, y0 = 0.0, ψ0 = 0.0)
    dim = length(time_vec)
    x_vec = zeros(Float64, dim)
    y_vec = zeros(Float64, dim)
    ψ_vec = zeros(Float64, dim)
    x_vec[1] = x0
    y_vec[1] = y0
    ψ_vec[1] = ψ0
    for i = 2:dim
        Δt = time_vec[i] - time_vec[i-1]
        ψ_vec[i] = ψ_vec[i-1] + r_vec[i] * Δt
        x_vec[i] = x_vec[i-1] + (u_vec[i] * cos(ψ_vec[i]) - v_vec[i] * sin(ψ_vec[i])) * Δt
        y_vec[i] = y_vec[i-1] + (u_vec[i] * sin(ψ_vec[i]) + v_vec[i] * cos(ψ_vec[i])) * Δt
    end
    x_vec, y_vec, ψ_vec
end

function estimate_kt_lsm(data::ShipData)
    time_vec = data.time
    r_vec = data.r
    δ_vec = data.δ

    spl_r = Spline1D(time_vec, r_vec)
    f = (spl, x) -> derivative(spl, x)
    dr_vec = f(spl_r, time_vec)
    estimate_kt_lsm(r_vec, dr_vec, δ_vec)
end

function estimate_kt_lsm(r_vec, dr_vec, δ_vec)
    A = hcat(δ_vec, r_vec)
    Θ = A \ dr_vec
    T = -1.0 / Θ[2]
    K = Θ[1] * T
    K, T
end

function estimate_kt_lsm_time_window_sampling(data::ShipData, window_size::Int)
    time_vec = data.time
    r = data.r
    spl_r = Spline1D(time_vec, r)
    f = (spl, x) -> derivative(spl, x)
    dr = f(spl_r, time_vec)
    δ = data.δ
    n_samples = length(time_vec) - window_size
    K_samples = zeros(n_samples)
    T_samples = zeros(n_samples)
    for i = 1:n_samples
        K_samples[i], T_samples[i] =
            estimate_kt_lsm(r[i:i+window_size], dr[i:i+window_size], δ[i:i+window_size])
    end
    K_samples, T_samples
end

function estimate_mmg_approx_lsm(
    data::ShipData,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2;
    ρ = 1.025,
)
    @unpack time, u, v, r, x, y, ψ, δ, npm = data
    @unpack L_pp,
    B,
    d,
    x_G,
    D_p,
    m,
    I_zG,
    A_R,
    η,
    m_x,
    m_y,
    J_z,
    f_α,
    ϵ,
    t_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    # -------------------------------------
    # calculating basic info
    U = sqrt.(u .^ 2 + (v - r .* x_G) .^ 2)

    β = U .^ 0
    for i = 1:length(time)
        if isapprox(U[i], 0.0) == false
            β[i] = asin(-(v[i] - r[i] * x_G) / U[i])
        end
    end

    v_dash = U .^ 0
    for i = 1:length(time)
        if isapprox(U[i], 0.0) == false
            v_dash[i] = v[i] / U[i]
        end
    end

    r_dash = U .^ 0
    for i = 1:length(time)
        if isapprox(U[i], 0.0) == false
            r_dash[i] = r[i] * L_pp / U[i]
        end
    end

    w_P = w_P0 * exp.(.-4.0 * (β - r_dash .* x_P) .^ 2)

    J = U .^ 0
    for i = 1:length(time)
        if isapprox(npm[i], 0.0) == false
            J[i] = (1.0 - w_P[i]) * u[i] / (npm[i] * D_p)
        end
    end

    K_T = J .^ 2 .* k_2 + J .* k_1 .+ k_0
    β_R = β - r_dash .* l_R

    γ_R = U .^ 0
    for i = 1:length(time)
        if β_R[i] < 0.0
            γ_R[i] = γ_R_minus
        else
            γ_R[i] = γ_R_plus
        end
    end

    v_R = U .* γ_R .* β_R

    u_R = U .^ 0
    for i = 1:length(time)
        if J[i] == 0.0
            u_R[i] = sqrt(η * (κ * ϵ * 8.0 * k_0 * npm[i]^2 * D_p^4 / pi)^2)
        else
            u_R[i] =
                u[i] *
                (1.0 - w_P[i]) *
                ϵ *
                sqrt(
                    η * (1.0 + κ * (sqrt(1.0 + 8.0 * K_T[i] / (pi * J[i]^2)) - 1.0))^2 +
                    (1.0 - η),
                )
        end
    end

    U_R = sqrt.(u_R .^ 2 + v_R .^ 2)
    α_R = δ - atan.(v_R, u_R)
    F_N = 0.5 * A_R * ρ * f_α * (U_R .^ 2) .* sin.(α_R)
    X_R = -(1.0 - t_R) .* F_N .* sin.(δ)
    X_P = (1.0 - t_P) * ρ .* K_T .* npm .^ 2 .* D_p^4
    Y_R = -(1.0 + a_H) .* F_N .* cos.(δ)
    N_R = -(-0.5 + a_H * x_H) .* F_N .* cos.(δ)
    # -------------------------------------

    # -------------------------------------
    spl_u = Spline1D(time, u)
    spl_v = Spline1D(time, v)
    spl_r = Spline1D(time, r)
    f = (spl, x) -> derivative(spl, x)
    du = f(spl_u, time)
    dv = f(spl_v, time)
    dr = f(spl_r, time)


    # calculating parameters
    Ax = hcat(
        -v_dash .^ 0,
        v_dash .^ 2,
        v_dash .* r_dash,
        r_dash .^ 2,
        (v_dash .^ 2) .* (r_dash .^ 2),
    )
    bx =
        ((m + m_x) .* du - (m + m_y) .* v .* r - x_G * m .* r .^ 2 - X_R - X_P) ./
        (0.5 * ρ * L_pp * d .* (U .^ 2))
    R_0_dash, X_vv_dash, X_vr_dash, X_rr_dash, X_vvvv_dash = Ax \ bx

    Ay = hcat(
        v_dash,
        r_dash,
        v_dash .^ 3,
        (v_dash .^ 2) .* (r_dash),
        (v_dash) .* (r_dash .^ 2),
        (r_dash .^ 3),
    )
    by =
        ((m + m_y) .* dv - (m + m_x) .* u .* r + x_G * m .* dr - Y_R) ./
        (0.5 * ρ * L_pp * d .* (U .^ 2))
    Y_v_dash, Y_r_dash, Y_vvv_dash, Y_vvr_dash, Y_vrr_dash, Y_rrr_dash = Ay \ by

    An = hcat(
        v_dash,
        r_dash,
        v_dash .^ 3,
        (v_dash .^ 2) .* (r_dash),
        (v_dash) .* (r_dash .^ 2),
        (r_dash .^ 3),
    )
    bn =
        ((I_zG + x_G^2 * m + J_z) .* dr + x_G * m .* (dv .+ u .* r) - N_R) ./
        (0.5 * ρ * L_pp .^ 2 * d .* (U .^ 2))
    N_v_dash, N_r_dash, N_vvv_dash, N_vvr_dash, N_vrr_dash, N_rrr_dash = An \ bn

    R_0_dash,
    X_vv_dash,
    X_vr_dash,
    X_rr_dash,
    X_vvvv_dash,
    Y_v_dash,
    Y_r_dash,
    Y_vvv_dash,
    Y_vvr_dash,
    Y_vrr_dash,
    Y_rrr_dash,
    N_v_dash,
    N_r_dash,
    N_vvv_dash,
    N_vvr_dash,
    N_vrr_dash,
    N_rrr_dash
end

function estimate_mmg_approx_lsm_time_window_sampling(
    data::ShipData,
    window_size::Int,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2;
    ρ = 1.025,
)
    n_samples = length(data.time) - window_size
    R_0_dash_samples = zeros(n_samples)
    X_vv_dash_samples = zeros(n_samples)
    X_vr_dash_samples = zeros(n_samples)
    X_rr_dash_samples = zeros(n_samples)
    X_vvvv_dash_samples = zeros(n_samples)
    Y_v_dash_samples = zeros(n_samples)
    Y_r_dash_samples = zeros(n_samples)
    Y_vvv_dash_samples = zeros(n_samples)
    Y_vvr_dash_samples = zeros(n_samples)
    Y_vrr_dash_samples = zeros(n_samples)
    Y_rrr_dash_samples = zeros(n_samples)
    N_v_dash_samples = zeros(n_samples)
    N_r_dash_samples = zeros(n_samples)
    N_vvv_dash_samples = zeros(n_samples)
    N_vvr_dash_samples = zeros(n_samples)
    N_vrr_dash_samples = zeros(n_samples)
    N_rrr_dash_samples = zeros(n_samples)
    for i = 1:n_samples
        sample_data = ShipData(
            data.time[i:i+window_size],
            data.u[i:i+window_size],
            data.v[i:i+window_size],
            data.r[i:i+window_size],
            data.x[i:i+window_size],
            data.y[i:i+window_size],
            data.ψ[i:i+window_size],
            data.δ[i:i+window_size],
            data.npm[i:i+window_size],
        )
        R_0_dash_samples[i],
        X_vv_dash_samples[i],
        X_vr_dash_samples[i],
        X_rr_dash_samples[i],
        X_vvvv_dash_samples[i],
        Y_v_dash_samples[i],
        Y_r_dash_samples[i],
        Y_vvv_dash_samples[i],
        Y_vvr_dash_samples[i],
        Y_vrr_dash_samples[i],
        Y_rrr_dash_samples[i],
        N_v_dash_samples[i],
        N_r_dash_samples[i],
        N_vvv_dash_samples[i],
        N_vvr_dash_samples[i],
        N_vrr_dash_samples[i],
        N_rrr_dash_samples[i] =
            estimate_mmg_approx_lsm(sample_data, basic_params, k_0, k_1, k_2, ρ = ρ)
    end
    R_0_dash_samples,
    X_vv_dash_samples,
    X_vr_dash_samples,
    X_rr_dash_samples,
    X_vvvv_dash_samples,
    Y_v_dash_samples,
    Y_r_dash_samples,
    Y_vvv_dash_samples,
    Y_vvr_dash_samples,
    Y_vrr_dash_samples,
    Y_rrr_dash_samples,
    N_v_dash_samples,
    N_r_dash_samples,
    N_vvv_dash_samples,
    N_vvr_dash_samples,
    N_vrr_dash_samples,
    N_rrr_dash_samples
end

function create_model_for_mcmc_sample_kt(
    data::ShipData;
    σ_r_prior_dist::Distribution = Chi(5),
    K_prior_dist::Distribution = Uniform(0.01, 10.0),
    T_prior_dist::Distribution = truncated(Normal(100.0, 50.0), 10.0, 200.0),
)
    time_obs = data.time
    r_obs = data.r
    δ_obs = data.δ

    # create sytem model
    spl_δ = Spline1D(time_obs, δ_obs)
    function KT!(dX, X, p, t)
        r, δ = X
        K, T = p
        dX[1] = dr = 1.0 / T * (-r + K * δ) # dr = 
        dX[2] = dδ = derivative(spl_δ, t) # dδ = 
    end
    u0 = [r_obs[1]; δ_obs[1]]
    K_start = 0.10
    T_start = 60.0
    p = [K_start, T_start]
    prob1 = ODEProblem(KT!, u0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    @model function fitKT(time_obs, r_obs, prob1)
        σ_r ~ σ_r_prior_dist
        K ~ K_prior_dist
        T ~ T_prior_dist

        p = [K, T]
        prob = remake(prob1, p = p)
        sol = solve(prob, Tsit5())
        predicted = sol(time_obs)
        for i = 1:length(predicted)
            r_obs[i] ~ Normal(predicted[i][1], σ_r) # index number of r is 1
        end
    end

    return fitKT(time_obs, r_obs, prob1)
end

function create_model_for_mcmc_sample_kt(
    data::ShipData,
    σ_r::Float64;
    K_prior_dist::Distribution = Uniform(0.01, 10.0),
    T_prior_dist::Distribution = truncated(Normal(100.0, 50.0), 10.0, 200.0),
)
    time_obs = data.time
    r_obs = data.r
    δ_obs = data.δ

    # create sytem model
    spl_δ = Spline1D(time_obs, δ_obs)
    function KT!(dX, X, p, t)
        r, δ = X
        K, T = p
        dX[1] = dr = 1.0 / T * (-r + K * δ) # dr = 
        dX[2] = dδ = derivative(spl_δ, t) # dδ = 
    end
    u0 = [r_obs[1]; δ_obs[1]]
    K_start = 0.10
    T_start = 60.0
    p = [K_start, T_start]
    prob1 = ODEProblem(KT!, u0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    @model function fitKT(time_obs, r_obs, prob1)
        K ~ K_prior_dist
        T ~ T_prior_dist

        p = [K, T]
        prob = remake(prob1, p = p)
        sol = solve(prob, Tsit5())
        predicted = sol(time_obs)
        for i = 1:length(predicted)
            r_obs[i] ~ Normal(predicted[i][1], σ_r) # index number of r is 1
        end
    end

    return fitKT(time_obs, r_obs, prob1)
end

function create_model_for_mcmc_sample_mmg(
    data::ShipData,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2;
    ρ = 1.025,
    σ_u_prior_dist = Uniform(0.00, 0.20),
    σ_v_prior_dist = Uniform(0.00, 0.20),
    σ_r_prior_dist = Uniform(0.00, 0.20),
    R_0_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vv_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vr_dash_prior_dist = Uniform(-2.0, 2.0),
    X_rr_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vvvv_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_v_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_r_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vvv_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vvr_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vrr_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_rrr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_v_dash_prior_dist = Uniform(-2.0, 2.0),
    N_r_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vvv_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vvr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vrr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_rrr_dash_prior_dist = Uniform(-2.0, 2.0),
    solver = Tsit5(),
    abstol = 1e-6,
    reltol = 1e-3,
)
    time_obs = data.time
    u_obs = data.u
    v_obs = data.v
    r_obs = data.r
    δ_obs = data.δ
    npm_obs = data.npm

    @unpack L_pp,
    B,
    d,
    x_G,
    D_p,
    m,
    I_zG,
    A_R,
    η,
    m_x,
    m_y,
    J_z,
    f_α,
    ϵ,
    t_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    # create sytem model
    spl_δ = Spline1D(time_obs, δ_obs)
    spl_npm = Spline1D(time_obs, npm_obs)
    function MMG!(dX, X, p, t)
        u, v, r, δ, npm = X
        R_0_dash,
        X_vv_dash,
        X_vr_dash,
        X_rr_dash,
        X_vvvv_dash,
        Y_v_dash,
        Y_r_dash,
        Y_vvv_dash,
        Y_vvr_dash,
        Y_vrr_dash,
        Y_rrr_dash,
        N_v_dash,
        N_r_dash,
        N_vvv_dash,
        N_vvr_dash,
        N_vrr_dash,
        N_rrr_dash = p

        U = sqrt(u^2 + (v - r * x_G)^2)

        β = 0.0
        if U == 0.0
            β = 0.0
        else
            β = asin(-(v - r * x_G) / U)
        end

        v_dash = 0.0
        if U == 0.0
            v_dash = 0.0
        else
            v_dash = v / U
        end

        r_dash = 0.0
        if U == 0.0
            r_dash = 0.0
        else
            r_dash = r * L_pp / U
        end

        w_P = w_P0 * exp(-4.0 * (β - x_P * r_dash)^2)

        J = 0.0
        if npm == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (npm * D_p)
        end
        K_T = k_0 + k_1 * J + k_2 * J^2
        β_R = β - l_R * r_dash
        γ_R = γ_R_minus

        if β_R < 0.0
            γ_R = γ_R_minus
        else
            γ_R = γ_R_plus
        end

        v_R = U * γ_R * β_R

        u_R = 0.0
        if J == 0.0
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * npm^2 * D_p^4 / pi)^2)
        else
            u_R =
                u *
                (1.0 - w_P) *
                ϵ *
                sqrt(η * (1.0 + κ * (sqrt(1.0 + 8.0 * K_T / (pi * J^2)) - 1))^2 + (1 - η))
        end

        U_R = sqrt(u_R^2 + v_R^2)
        α_R = δ - atan(v_R, u_R)
        F_N = 0.5 * A_R * ρ * f_α * (U_R^2) * sin(α_R)

        X_H = (
            0.5 *
            ρ *
            L_pp *
            d *
            (U^2) *
            (
                -R_0_dash +
                X_vv_dash * (v_dash^2) +
                X_vr_dash * v_dash * r_dash +
                X_rr_dash * (r_dash^2) +
                X_vvvv_dash * (v_dash^4)
            )
        )
        X_R = -(1.0 - t_R) * F_N * sin(δ)
        X_P = (1.0 - t_P) * ρ * K_T * npm^2 * D_p^4
        Y_H = (
            0.5 *
            ρ *
            L_pp *
            d *
            (U^2) *
            (
                Y_v_dash * v_dash +
                Y_r_dash * r_dash +
                Y_vvv_dash * (v_dash^3) +
                Y_vvr_dash * (v_dash^2) * r_dash +
                Y_vrr_dash * v_dash * (r_dash^2) +
                Y_rrr_dash * (r_dash^3)
            )
        )
        Y_R = -(1 + a_H) * F_N * cos(δ)
        N_H = (
            0.5 *
            ρ *
            (L_pp^2) *
            d *
            (U^2) *
            (
                N_v_dash * v_dash +
                N_r_dash * r_dash +
                N_vvv_dash * (v_dash^3) +
                N_vvr_dash * (v_dash^2) * r_dash +
                N_vrr_dash * v_dash * (r_dash^2) +
                N_rrr_dash * (r_dash^3)
            )
        )
        N_R = -(-0.5 + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                    ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dδ = derivative(spl_δ, t)
        dX[5] = dnpm = derivative(spl_npm, t)
    end

    R_0_dash_start = 0.022
    X_vv_dash_start = -0.040
    X_vr_dash_start = 0.002
    X_rr_dash_start = 0.011
    X_vvvv_dash_start = 0.771
    Y_v_dash_start = -0.315
    Y_r_dash_start = 0.083
    Y_vvv_dash_start = -1.607
    Y_vvr_dash_start = 0.379
    Y_vrr_dash_start = -0.391
    Y_rrr_dash_start = 0.008
    N_v_dash_start = -0.137
    N_r_dash_start = -0.049
    N_vvv_dash_start = -0.030
    N_vvr_dash_start = -0.294
    N_vrr_dash_start = 0.055
    N_rrr_dash_start = -0.013

    p = [
        R_0_dash_start,
        X_vv_dash_start,
        X_vr_dash_start,
        X_rr_dash_start,
        X_vvvv_dash_start,
        Y_v_dash_start,
        Y_r_dash_start,
        Y_vvv_dash_start,
        Y_vvr_dash_start,
        Y_vrr_dash_start,
        Y_rrr_dash_start,
        N_v_dash_start,
        N_r_dash_start,
        N_vvv_dash_start,
        N_vvr_dash_start,
        N_vrr_dash_start,
        N_rrr_dash_start,
    ]

    u0 = 2.29 * 0.512
    v0 = 0.0
    r0 = 0.0
    X0 = [u_obs[1]; v_obs[1]; r_obs[1]; δ_obs[1]; npm_obs[1]]
    prob1 = ODEProblem(MMG!, X0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    @model function fitMMG(time_obs, obs, prob1)
        σ_u ~ σ_u_prior_dist
        σ_v ~ σ_v_prior_dist
        σ_r ~ σ_r_prior_dist
        R_0_dash ~ R_0_dash_prior_dist
        X_vv_dash ~ X_vv_dash_prior_dist
        X_vr_dash ~ X_vr_dash_prior_dist
        X_rr_dash ~ X_rr_dash_prior_dist
        X_vvvv_dash ~ X_vvvv_dash_prior_dist
        Y_v_dash ~ Y_v_dash_prior_dist
        Y_r_dash ~ Y_r_dash_prior_dist
        Y_vvv_dash ~ Y_vvv_dash_prior_dist
        Y_vvr_dash ~ Y_vvr_dash_prior_dist
        Y_vrr_dash ~ Y_vrr_dash_prior_dist
        Y_rrr_dash ~ Y_rrr_dash_prior_dist
        N_v_dash ~ N_v_dash_prior_dist
        N_r_dash ~ N_r_dash_prior_dist
        N_vvv_dash ~ N_vvv_dash_prior_dist
        N_vvr_dash ~ N_vvr_dash_prior_dist
        N_vrr_dash ~ N_vrr_dash_prior_dist
        N_rrr_dash ~ N_rrr_dash_prior_dist

        p = [
            R_0_dash,
            X_vv_dash,
            X_vr_dash,
            X_rr_dash,
            X_vvvv_dash,
            Y_v_dash,
            Y_r_dash,
            Y_vvv_dash,
            Y_vvr_dash,
            Y_vrr_dash,
            Y_rrr_dash,
            N_v_dash,
            N_r_dash,
            N_vvv_dash,
            N_vvr_dash,
            N_vrr_dash,
            N_rrr_dash,
        ]
        prob = remake(prob1, p = p)
        sol = solve(prob, solver, abstol = abstol, reltol = reltol)
        predicted = sol(time_obs)
        for i = 1:length(predicted)
            obs[1][i] ~ Normal(predicted[i][1], σ_u) # u
            obs[2][i] ~ Normal(predicted[i][2], σ_v) # v
            obs[3][i] ~ Normal(predicted[i][3], σ_r) # r
        end
    end

    return fitMMG(time_obs, [u_obs, v_obs, r_obs], prob1)
end

function create_model_for_mcmc_sample_mmg(
    data::ShipData,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2,
    σ_u,
    σ_v,
    σ_r;
    ρ = 1.025,
    R_0_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vv_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vr_dash_prior_dist = Uniform(-2.0, 2.0),
    X_rr_dash_prior_dist = Uniform(-2.0, 2.0),
    X_vvvv_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_v_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_r_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vvv_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vvr_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_vrr_dash_prior_dist = Uniform(-2.0, 2.0),
    Y_rrr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_v_dash_prior_dist = Uniform(-2.0, 2.0),
    N_r_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vvv_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vvr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_vrr_dash_prior_dist = Uniform(-2.0, 2.0),
    N_rrr_dash_prior_dist = Uniform(-2.0, 2.0),
    solver = Tsit5(),
    abstol = 1e-6,
    reltol = 1e-3,
)
    time_obs = data.time
    u_obs = data.u
    v_obs = data.v
    r_obs = data.r
    δ_obs = data.δ
    npm_obs = data.npm

    @unpack L_pp,
    B,
    d,
    x_G,
    D_p,
    m,
    I_zG,
    A_R,
    η,
    m_x,
    m_y,
    J_z,
    f_α,
    ϵ,
    t_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    # create sytem model
    spl_δ = Spline1D(time_obs, δ_obs)
    spl_npm = Spline1D(time_obs, npm_obs)
    function MMG!(dX, X, p, t)
        u, v, r, δ, npm = X
        R_0_dash,
        X_vv_dash,
        X_vr_dash,
        X_rr_dash,
        X_vvvv_dash,
        Y_v_dash,
        Y_r_dash,
        Y_vvv_dash,
        Y_vvr_dash,
        Y_vrr_dash,
        Y_rrr_dash,
        N_v_dash,
        N_r_dash,
        N_vvv_dash,
        N_vvr_dash,
        N_vrr_dash,
        N_rrr_dash = p

        U = sqrt(u^2 + (v - r * x_G)^2)

        β = 0.0
        if U == 0.0
            β = 0.0
        else
            β = asin(-(v - r * x_G) / U)
        end

        v_dash = 0.0
        if U == 0.0
            v_dash = 0.0
        else
            v_dash = v / U
        end

        r_dash = 0.0
        if U == 0.0
            r_dash = 0.0
        else
            r_dash = r * L_pp / U
        end

        w_P = w_P0 * exp(-4.0 * (β - x_P * r_dash)^2)

        J = 0.0
        if npm == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (npm * D_p)
        end
        K_T = k_0 + k_1 * J + k_2 * J^2
        β_R = β - l_R * r_dash
        γ_R = γ_R_minus

        if β_R < 0.0
            γ_R = γ_R_minus
        else
            γ_R = γ_R_plus
        end

        v_R = U * γ_R * β_R

        u_R = 0.0
        if J == 0.0
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * npm^2 * D_p^4 / pi)^2)
        else
            u_R =
                u *
                (1.0 - w_P) *
                ϵ *
                sqrt(η * (1.0 + κ * (sqrt(1.0 + 8.0 * K_T / (pi * J^2)) - 1))^2 + (1 - η))
        end

        U_R = sqrt(u_R^2 + v_R^2)
        α_R = δ - atan(v_R, u_R)
        F_N = 0.5 * A_R * ρ * f_α * (U_R^2) * sin(α_R)

        X_H = (
            0.5 *
            ρ *
            L_pp *
            d *
            (U^2) *
            (
                -R_0_dash +
                X_vv_dash * (v_dash^2) +
                X_vr_dash * v_dash * r_dash +
                X_rr_dash * (r_dash^2) +
                X_vvvv_dash * (v_dash^4)
            )
        )
        X_R = -(1.0 - t_R) * F_N * sin(δ)
        X_P = (1.0 - t_P) * ρ * K_T * npm^2 * D_p^4
        Y_H = (
            0.5 *
            ρ *
            L_pp *
            d *
            (U^2) *
            (
                Y_v_dash * v_dash +
                Y_r_dash * r_dash +
                Y_vvv_dash * (v_dash^3) +
                Y_vvr_dash * (v_dash^2) * r_dash +
                Y_vrr_dash * v_dash * (r_dash^2) +
                Y_rrr_dash * (r_dash^3)
            )
        )
        Y_R = -(1 + a_H) * F_N * cos(δ)
        N_H = (
            0.5 *
            ρ *
            (L_pp^2) *
            d *
            (U^2) *
            (
                N_v_dash * v_dash +
                N_r_dash * r_dash +
                N_vvv_dash * (v_dash^3) +
                N_vvr_dash * (v_dash^2) * r_dash +
                N_vrr_dash * v_dash * (r_dash^2) +
                N_rrr_dash * (r_dash^3)
            )
        )
        N_R = -(-0.5 + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                    ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dδ = derivative(spl_δ, t)
        dX[5] = dnpm = derivative(spl_npm, t)
    end

    R_0_dash_start = 0.022
    X_vv_dash_start = -0.040
    X_vr_dash_start = 0.002
    X_rr_dash_start = 0.011
    X_vvvv_dash_start = 0.771
    Y_v_dash_start = -0.315
    Y_r_dash_start = 0.083
    Y_vvv_dash_start = -1.607
    Y_vvr_dash_start = 0.379
    Y_vrr_dash_start = -0.391
    Y_rrr_dash_start = 0.008
    N_v_dash_start = -0.137
    N_r_dash_start = -0.049
    N_vvv_dash_start = -0.030
    N_vvr_dash_start = -0.294
    N_vrr_dash_start = 0.055
    N_rrr_dash_start = -0.013

    p = [
        R_0_dash_start,
        X_vv_dash_start,
        X_vr_dash_start,
        X_rr_dash_start,
        X_vvvv_dash_start,
        Y_v_dash_start,
        Y_r_dash_start,
        Y_vvv_dash_start,
        Y_vvr_dash_start,
        Y_vrr_dash_start,
        Y_rrr_dash_start,
        N_v_dash_start,
        N_r_dash_start,
        N_vvv_dash_start,
        N_vvr_dash_start,
        N_vrr_dash_start,
        N_rrr_dash_start,
    ]

    u0 = 2.29 * 0.512
    v0 = 0.0
    r0 = 0.0
    X0 = [u_obs[1]; v_obs[1]; r_obs[1]; δ_obs[1]; npm_obs[1]]
    prob1 = ODEProblem(MMG!, X0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    @model function fitMMG(time_obs, obs, prob1)
        R_0_dash ~ R_0_dash_prior_dist
        X_vv_dash ~ X_vv_dash_prior_dist
        X_vr_dash ~ X_vr_dash_prior_dist
        X_rr_dash ~ X_rr_dash_prior_dist
        X_vvvv_dash ~ X_vvvv_dash_prior_dist
        Y_v_dash ~ Y_v_dash_prior_dist
        Y_r_dash ~ Y_r_dash_prior_dist
        Y_vvv_dash ~ Y_vvv_dash_prior_dist
        Y_vvr_dash ~ Y_vvr_dash_prior_dist
        Y_vrr_dash ~ Y_vrr_dash_prior_dist
        Y_rrr_dash ~ Y_rrr_dash_prior_dist
        N_v_dash ~ N_v_dash_prior_dist
        N_r_dash ~ N_r_dash_prior_dist
        N_vvv_dash ~ N_vvv_dash_prior_dist
        N_vvr_dash ~ N_vvr_dash_prior_dist
        N_vrr_dash ~ N_vrr_dash_prior_dist
        N_rrr_dash ~ N_rrr_dash_prior_dist

        p = [
            R_0_dash,
            X_vv_dash,
            X_vr_dash,
            X_rr_dash,
            X_vvvv_dash,
            Y_v_dash,
            Y_r_dash,
            Y_vvv_dash,
            Y_vvr_dash,
            Y_vrr_dash,
            Y_rrr_dash,
            N_v_dash,
            N_r_dash,
            N_vvv_dash,
            N_vvr_dash,
            N_vrr_dash,
            N_rrr_dash,
        ]
        prob = remake(prob1, p = p)
        sol = solve(prob, solver, abstol = abstol, reltol = reltol)
        predicted = sol(time_obs)
        for i = 1:length(predicted)
            obs[1][i] ~ Normal(predicted[i][1], σ_u) # u
            obs[2][i] ~ Normal(predicted[i][2], σ_v) # v
            obs[3][i] ~ Normal(predicted[i][3], σ_r) # r
        end
    end

    return fitMMG(time_obs, [u_obs, v_obs, r_obs], prob1)
end

function nuts_sampling_single_thread(
    model,
    n_samples::Int,
    n_chains::Int;
    target_acceptance::Float64 = 0.65,
    progress = true,
)
    sampler = NUTS(target_acceptance)
    mapreduce(
        c -> sample(model, sampler, n_samples, progress = progress),
        chainscat,
        1:n_chains,
    )
end

function nuts_sampling_multi_threads(
    model,
    n_samples::Int,
    n_chains::Int;
    target_acceptance::Float64 = 0.65,
    progress = false,
)
    sampler = NUTS(target_acceptance)
    sample(model, sampler, MCMCThreads(), n_samples, n_chains, progress = progress)
end