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

@testset "kt bootstrap least square method func" begin
    K = 0.155  # [1/s]
    T = 80.5  # [s]
    duration = 500  # [s]
    sampling = 10000
    time_list = range(0.0, stop = duration, length = sampling)
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    r, δ = kt_simulate(K, T, time_list, δ_list)
    noize_dist = Normal(0.0, 0.0005)
    r_obs = r + rand(noize_dist, size(r))
    data = ShipData(time_list, 0, 0, r_obs, 0, 0, 0, δ, 0)
    one_sample_size = 5000
    K_est_samples, T_est_samples =
        estimate_kt_lsm_time_window_sampling(data, one_sample_size)
end

@testset "kt mcmc sampling func" begin
    K = 0.155  # [1/s]
    T = 80.5  # [s]
    duration = 500  # [s]
    sampling = 1001
    time_list = range(0.0, stop = duration, length = sampling)
    Ts = 50.0
    δ_list = 10.0 * pi / 180.0 * sin.(2.0 * pi / Ts * time_list) # [rad]
    r, δ = kt_simulate(K, T, time_list, δ_list)
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

@testset "mmg 3dof approx least square method func" begin
    # --- KVLCC2_L7 --
    ρ = 1.025  # 海水密度

    L_pp = 7.00  # 船長Lpp[m]
    B = 1.27  # 船幅[m]
    d = 0.46  # 喫水[m]
    nabla = 3.27  # 排水量[m^3]
    x_G = 0.25  # 重心位置[m]
    # C_b = 0.810  # 方形係数[-]
    D_p = 0.216  # プロペラ直径[m]
    H_R = 0.345  # 舵高さ[m]
    A_R = 0.0539  # 舵断面積[m^2]
    t_P = 0.220  # 推力減少率
    w_P0 = 0.40  # 有効伴流率
    m_x_dash = 0.022  # 付加質量x(無次元)
    m_y_dash = 0.223  # 付加質量y(無次元)
    J_z_dash = 0.011  # 付加質量Izz(無次元)
    t_R = 0.387  # 操縦抵抗減少率
    a_H = 0.312  # 舵力増加係数
    x_H_dash = -0.464  # 舵力増分作用位置
    γ_R_minus = 0.395  # 整流係数
    γ_R_plus = 0.640  # 整流係数
    l_r_dash = -0.710  # 船長に対する舵位置
    x_P_dash = -0.690  # 船長に対するプロペラ位置
    ϵ = 1.09  # プロペラ・舵位置伴流係数比
    κ = 0.50  # 修正係数
    f_α = 2.747  # 直圧力勾配係数

    L_pp = L_pp  # 船長Lpp[m]
    B = B  # 船幅[m]
    d = d  # 喫水[m]
    x_G = x_G  # 重心位置[]
    D_p = D_p  # プロペラ直径[m]
    m = ρ * nabla  # 質量(無次元化)[kg]
    I_zG = ρ * nabla * ((0.25 * L_pp)^2)  # 慣性モーメント[-]
    A_R = A_R  # 船の断面に対する舵面積比[-]
    η = D_p / H_R  # プロペラ直径に対する舵高さ(Dp/H)
    m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash  # 付加質量x(無次元)
    m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash  # 付加質量y(無次元)
    J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash  # 付加質量Izz(無次元)
    f_α = f_α # 直圧力勾配係数
    ϵ = ϵ  # プロペラ・舵位置伴流係数比
    t_R = t_R  # 操縦抵抗減少率
    a_H = a_H  # 舵力増加係数
    x_H = x_H_dash * L_pp  # 舵力増分作用位置
    γ_R_minus = γ_R_minus  # 整流係数
    γ_R_plus = γ_R_plus  # 整流係数
    l_R = l_r_dash  # 船長に対する舵位置
    κ = κ  # 修正係数
    t_P = t_P  # 推力減少率
    w_P0 = w_P0  # 有効伴流率
    x_P = x_P_dash  # 船長に対するプロペラ位置
    basic_params = Mmg3DofBasicParams(
        L_pp,
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
        x_P,
    )

    k_0 = 0.2931
    k_1 = -0.2753
    k_2 = -0.1385
    R_0_dash = 0.022
    X_vv_dash = -0.040
    X_vr_dash = 0.002
    X_rr_dash = 0.011
    X_vvvv_dash = 0.771
    Y_v_dash = -0.315
    Y_r_dash = 0.083
    Y_vvv_dash = -1.607
    Y_vvr_dash = 0.379
    Y_vrr_dash = -0.391
    Y_rrr_dash = 0.008
    N_v_dash = -0.137
    N_r_dash = -0.049
    N_vvv_dash = -0.030
    N_vvr_dash = -0.294
    N_vrr_dash = 0.055
    N_rrr_dash = -0.013
    maneuvering_params = Mmg3DofManeuveringParams(
        k_0,
        k_1,
        k_2,
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
    )
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rpm]

    sampling = duration * 10
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    npm_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        npm_list,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    u, v, r, δ, npm = mmg_results
    x, y, ψ = calc_position(time_list, u, v, r)
    data = ShipData(time_list, u, v, r, x, y, ψ, δ, npm)
    estimate_mmg_approx_lsm(data)
end

@testset "mmg 3dof approx time window sampling least square method func" begin
    # --- KVLCC2_L7 --
    ρ = 1.025  # 海水密度

    L_pp = 7.00  # 船長Lpp[m]
    B = 1.27  # 船幅[m]
    d = 0.46  # 喫水[m]
    nabla = 3.27  # 排水量[m^3]
    x_G = 0.25  # 重心位置[m]
    # C_b = 0.810  # 方形係数[-]
    D_p = 0.216  # プロペラ直径[m]
    H_R = 0.345  # 舵高さ[m]
    A_R = 0.0539  # 舵断面積[m^2]
    t_P = 0.220  # 推力減少率
    w_P0 = 0.40  # 有効伴流率
    m_x_dash = 0.022  # 付加質量x(無次元)
    m_y_dash = 0.223  # 付加質量y(無次元)
    J_z_dash = 0.011  # 付加質量Izz(無次元)
    t_R = 0.387  # 操縦抵抗減少率
    a_H = 0.312  # 舵力増加係数
    x_H_dash = -0.464  # 舵力増分作用位置
    γ_R_minus = 0.395  # 整流係数
    γ_R_plus = 0.640  # 整流係数
    l_r_dash = -0.710  # 船長に対する舵位置
    x_P_dash = -0.690  # 船長に対するプロペラ位置
    ϵ = 1.09  # プロペラ・舵位置伴流係数比
    κ = 0.50  # 修正係数
    f_α = 2.747  # 直圧力勾配係数

    L_pp = L_pp  # 船長Lpp[m]
    B = B  # 船幅[m]
    d = d  # 喫水[m]
    x_G = x_G  # 重心位置[]
    D_p = D_p  # プロペラ直径[m]
    m = ρ * nabla  # 質量(無次元化)[kg]
    I_zG = ρ * nabla * ((0.25 * L_pp)^2)  # 慣性モーメント[-]
    A_R = A_R  # 船の断面に対する舵面積比[-]
    η = D_p / H_R  # プロペラ直径に対する舵高さ(Dp/H)
    m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash  # 付加質量x(無次元)
    m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash  # 付加質量y(無次元)
    J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash  # 付加質量Izz(無次元)
    f_α = f_α # 直圧力勾配係数
    ϵ = ϵ  # プロペラ・舵位置伴流係数比
    t_R = t_R  # 操縦抵抗減少率
    a_H = a_H  # 舵力増加係数
    x_H = x_H_dash * L_pp  # 舵力増分作用位置
    γ_R_minus = γ_R_minus  # 整流係数
    γ_R_plus = γ_R_plus  # 整流係数
    l_R = l_r_dash  # 船長に対する舵位置
    κ = κ  # 修正係数
    t_P = t_P  # 推力減少率
    w_P0 = w_P0  # 有効伴流率
    x_P = x_P_dash  # 船長に対するプロペラ位置
    basic_params = Mmg3DofBasicParams(
        L_pp,
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
        x_P,
    )

    k_0 = 0.2931
    k_1 = -0.2753
    k_2 = -0.1385
    R_0_dash = 0.022
    X_vv_dash = -0.040
    X_vr_dash = 0.002
    X_rr_dash = 0.011
    X_vvvv_dash = 0.771
    Y_v_dash = -0.315
    Y_r_dash = 0.083
    Y_vvv_dash = -1.607
    Y_vvr_dash = 0.379
    Y_vrr_dash = -0.391
    Y_rrr_dash = 0.008
    N_v_dash = -0.137
    N_r_dash = -0.049
    N_vvv_dash = -0.030
    N_vvr_dash = -0.294
    N_vrr_dash = 0.055
    N_rrr_dash = -0.013
    maneuvering_params = Mmg3DofManeuveringParams(
        k_0,
        k_1,
        k_2,
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
    )
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rpm]

    sampling = duration * 10
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    npm_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        npm_list,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    u, v, r, δ, npm = mmg_results
    noize_dist = Normal(0.0, 0.0005)
    u_obs = u + rand(noize_dist, size(u))
    v_obs = v + rand(noize_dist, size(v))
    r_obs = r + rand(noize_dist, size(r))
    x, y, ψ = calc_position(time_list, u_obs, v_obs, r_obs)
    data = ShipData(time_list, u_obs, v_obs, r_obs, x, y, ψ, δ, npm)
    one_sample_size = 5000
    estimate_mmg_apporx_lsm_time_window_sampling(data, 5000)
end

@testset "mmg mcmc sampling func" begin
    # --- KVLCC2_L7 --
    ρ = 1.025  # 海水密度

    L_pp = 7.00  # 船長Lpp[m]
    B = 1.27  # 船幅[m]
    d = 0.46  # 喫水[m]
    nabla = 3.27  # 排水量[m^3]
    x_G = 0.25  # 重心位置[m]
    # C_b = 0.810  # 方形係数[-]
    D_p = 0.216  # プロペラ直径[m]
    H_R = 0.345  # 舵高さ[m]
    A_R = 0.0539  # 舵断面積[m^2]
    t_P = 0.220  # 推力減少率
    w_P0 = 0.40  # 有効伴流率
    m_x_dash = 0.022  # 付加質量x(無次元)
    m_y_dash = 0.223  # 付加質量y(無次元)
    J_z_dash = 0.011  # 付加質量Izz(無次元)
    t_R = 0.387  # 操縦抵抗減少率
    a_H = 0.312  # 舵力増加係数
    x_H_dash = -0.464  # 舵力増分作用位置
    γ_R_minus = 0.395  # 整流係数
    γ_R_plus = 0.640  # 整流係数
    l_r_dash = -0.710  # 船長に対する舵位置
    x_P_dash = -0.690  # 船長に対するプロペラ位置
    ϵ = 1.09  # プロペラ・舵位置伴流係数比
    κ = 0.50  # 修正係数
    f_α = 2.747  # 直圧力勾配係数

    L_pp = L_pp  # 船長Lpp[m]
    B = B  # 船幅[m]
    d = d  # 喫水[m]
    x_G = x_G  # 重心位置[]
    D_p = D_p  # プロペラ直径[m]
    m = ρ * nabla  # 質量(無次元化)[kg]
    I_zG = ρ * nabla * ((0.25 * L_pp)^2)  # 慣性モーメント[-]
    A_R = A_R  # 船の断面に対する舵面積比[-]
    η = D_p / H_R  # プロペラ直径に対する舵高さ(Dp/H)
    m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash  # 付加質量x(無次元)
    m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash  # 付加質量y(無次元)
    J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash  # 付加質量Izz(無次元)
    f_α = f_α # 直圧力勾配係数
    ϵ = ϵ  # プロペラ・舵位置伴流係数比
    t_R = t_R  # 操縦抵抗減少率
    a_H = a_H  # 舵力増加係数
    x_H = x_H_dash * L_pp  # 舵力増分作用位置
    γ_R_minus = γ_R_minus  # 整流係数
    γ_R_plus = γ_R_plus  # 整流係数
    l_R = l_r_dash  # 船長に対する舵位置
    κ = κ  # 修正係数
    t_P = t_P  # 推力減少率
    w_P0 = w_P0  # 有効伴流率
    x_P = x_P_dash  # 船長に対するプロペラ位置
    basic_params = Mmg3DofBasicParams(
        L_pp,
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
        x_P,
    )

    k_0 = 0.2931
    k_1 = -0.2753
    k_2 = -0.1385
    R_0_dash = 0.022
    X_vv_dash = -0.040
    X_vr_dash = 0.002
    X_rr_dash = 0.011
    X_vvvv_dash = 0.771
    Y_v_dash = -0.315
    Y_r_dash = 0.083
    Y_vvv_dash = -1.607
    Y_vvr_dash = 0.379
    Y_vrr_dash = -0.391
    Y_rrr_dash = 0.008
    N_v_dash = -0.137
    N_r_dash = -0.049
    N_vvv_dash = -0.030
    N_vvr_dash = -0.294
    N_vrr_dash = 0.055
    N_rrr_dash = -0.013
    maneuvering_params = Mmg3DofManeuveringParams(
        k_0,
        k_1,
        k_2,
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
    )
    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rpm]

    sampling = duration * 10 + 1
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    npm_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        basic_params,
        maneuvering_params,
        time_list,
        δ_rad_list,
        npm_list,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    u, v, r, δ, npm = mmg_results
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
        npm[1:sampling_rate:end],
    )
    n_samples = 10
    n_chains = 1
    model = create_model_for_mcmc_sample_mmg(data, basic_params, k_0, k_1, k_2)
    chain = nuts_sampling_single_thread(model, n_samples, n_chains)
end
