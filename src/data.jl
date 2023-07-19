@with_kw struct ShipData{Tt,Tu,Tv,Tr,Tx,Ty,Tψ,Tδ,Tn_p}
    time::Tt
    u::Tu
    v::Tv
    r::Tr
    x::Tx
    y::Ty
    ψ::Tψ
    δ::Tδ
    n_p::Tn_p
end

function get_KVLCC2_L7_basic_params(ρ=1025.0)
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
    x_R_dash = -0.500  # 舵の相対位置
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
    x_R = x_R_dash * L_pp  # 舵の位置
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
        x_R,
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
    basic_params
end

function get_KVLCC2_L7_maneuvering_params()
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
    maneuvering_params
end

function get_KVLCC2_L7_params()
    basic_params = get_KVLCC2_L7_basic_params()
    maneuvering_params = get_KVLCC2_L7_maneuvering_params()
    basic_params, maneuvering_params
end

"""
    wind_force_and_moment_coefficients(ψ_A,p)

compute the parameter C_X,C_Y,C_N from Fujiwara's estimation method (1994).

# Arguments
-`ψ_A`: Wind attack of angle [rad]
- L_pp
- B
- A_OD
- A_F
- A_L
- H_BR
- H_C
- C
"""
function wind_force_and_moment_coefficients(
    ψ_A,
    L_pp,
    B,
    A_OD,
    A_F,
    A_L,
    H_BR,
    H_C,
    C,
)
    #C_LF1の場合で調整
    C_CF = 0.404 + 0.368 * A_F / (B * H_BR) + 0.902 * H_BR / L_pp

    if deg2rad(0) <= ψ_A <= deg2rad(90)
        C_LF = -0.992 + 0.507 * A_L / (L_pp * B) + 1.162 * C / L_pp
        C_XLI = 0.458 + 3.245 * A_L / (L_pp * H_BR) - 2.313 * A_F / (B * H_BR)
        C_ALF = -0.585 - 0.906 * A_OD / A_L + 3.239 * B / L_pp
        C_YLI = pi * A_L / L_pp^2 + 0.116 + 3.345 * A_F / (L_pp * B)

    elseif deg2rad(90) < ψ_A <= deg2rad(180)
        C_LF =
            0.018 - 5.091 * B / L_pp + 10.367 * H_C / L_pp - 3.011 * A_OD / L_pp^2 -
            0.341 * A_F / B^2
        C_XLI =
            -1.901 + 12.727 * A_L / (L_pp * H_BR) + 24.407 * A_F / A_L -
            40.310 * B / L_pp - 0.341 * A_F / (B * H_BR)
        C_ALF = -0.314 - 1.117 * A_OD / A_L
        C_YLI = pi * A_L / L_pp^2 + 0.446 + 2.192 * A_F / L_pp^2

    elseif deg2rad(180) < ψ_A <= deg2rad(270)
        C_LF =
            0.018 - 5.091 * B / L_pp + 10.367 * H_C / L_pp - 3.011 * A_OD / L_pp^2 -
            0.341 * A_F / B^2
        C_XLI =
            -1.901 + 12.727 * A_L / (L_pp * H_BR) + 24.407 * A_F / A_L -
            40.310 * B / L_pp - 0.341 * A_F / (B * H_BR)
        C_ALF = -(-0.314 - 1.117 * A_OD / A_L)
        C_YLI = -(pi * A_L / L_pp^2 + 0.446 + 2.192 * A_F / L_pp^2)

    elseif deg2rad(270) < ψ_A <= deg2rad(360)
        C_LF = -0.992 + 0.507 * A_L / (L_pp * B) + 1.162 * C / L_pp
        C_XLI = 0.458 + 3.245 * A_L / (L_pp * H_BR) - 2.313 * A_F / (B * H_BR)
        C_ALF = -(-0.585 - 0.906 * A_OD / A_L + 3.239 * B / L_pp)
        C_YLI = -(pi * A_L / L_pp^2 + 0.116 + 3.345 * A_F / (L_pp * B))
    end

    C_X =
        C_LF * cos(ψ_A) +
        C_XLI * (sin(ψ_A) - sin(ψ_A) * cos(ψ_A)^2 / 2) * sin(ψ_A) * cos(ψ_A) +
        C_ALF * sin(ψ_A) * cos(ψ_A)^3
    C_Y =
        C_CF * sin(ψ_A)^2 +
        C_YLI * (cos(ψ_A) + sin(ψ_A)^2 * cos(ψ_A) / 2) * sin(ψ_A) * cos(ψ_A)
    C_N = C_Y * (0.297 * C / L_pp - 0.149 * (ψ_A - deg2rad(90)))

    C_X, C_Y, C_N
end

function get_example_ship_wind_force_moment_params()
    L_pp = 7.00  # 船長Lpp[m]
    B = 1.27  # 船幅[m]
    d = 0.46  # 喫水[m]
    D = 0.6563 # 深さ[m]
    A_OD = 0.65 # デッキ上の構造物の側面投影面積[m^2]
    H_BR = 0.85 # 喫水からブリッジ主要構造物の最高位[m]
    H_C = 0.235 # 喫水から側面積中心までの高さ[m]
    C = 0.0 # 船体中心から側面積中心までの前後方向座標(船首方向を正)[m]

    A_OD = A_OD # デッキ上の構造物の側面投影面積[m^2]
    A_F = (D - d) * B  # 船体の正面投影面積[m^2]
    A_L = (D - d) * L_pp # 船体の側面投影面積[m^2]
    H_BR = H_BR # 喫水からブリッジ主要構造物の最高位[m]
    H_C = H_C # 喫水から側面積中心までの高さ[m]
    C = C # 船体中心から側面積中心までの前後方向座標[m]

    ψ_A_vec = deg2rad.(collect(0:10:360))
    C_X_vec = Array{Float64}(undef, length(ψ_A_vec))
    C_Y_vec = Array{Float64}(undef, length(ψ_A_vec))
    C_N_vec = Array{Float64}(undef, length(ψ_A_vec))
    for (index, ψ_A) in enumerate(ψ_A_vec)
        C_X, C_Y, C_N = wind_force_and_moment_coefficients(
            ψ_A,
            L_pp,
            B,
            A_OD,
            A_F,
            A_L,
            H_BR,
            H_C,
            C,
        )
        C_X_vec[index] = C_X
        C_Y_vec[index] = C_Y
        C_N_vec[index] = C_N
    end
    spl_C_X = Spline1D(ψ_A_vec, C_X_vec)
    spl_C_Y = Spline1D(ψ_A_vec, C_Y_vec)
    spl_C_N = Spline1D(ψ_A_vec, C_N_vec)

    Mmg3DofWindForceMomentParams(A_F, A_L, spl_C_X, spl_C_Y, spl_C_N)
end

function calc_position(time_vec, u_vec, v_vec, r_vec; x0=0.0, y0=0.0, ψ0=0.0)
    dim = length(time_vec)
    x_vec = zeros(Float64, dim)
    y_vec = zeros(Float64, dim)
    ψ_vec = zeros(Float64, dim)
    x_vec[1] = x0
    y_vec[1] = y0
    ψ_vec[1] = ψ0
    for i in 2:dim
        Δt = time_vec[i] - time_vec[i-1]
        ψ_vec[i] = ψ_vec[i-1] + r_vec[i] * Δt
        x_vec[i] =
            x_vec[i-1] + (u_vec[i] * cos(ψ_vec[i]) - v_vec[i] * sin(ψ_vec[i])) * Δt
        y_vec[i] =
            y_vec[i-1] + (u_vec[i] * sin(ψ_vec[i]) + v_vec[i] * cos(ψ_vec[i])) * Δt
    end
    x_vec, y_vec, ψ_vec
end

function nuts_sampling_single_thread(
    model,
    n_samples::Int,
    n_chains::Int;
    target_acceptance::Float64=0.65,
    progress=true,
)
    sampler = NUTS(target_acceptance)
    mapreduce(
        c -> sample(model, sampler, n_samples, progress=progress),
        chainscat,
        1:n_chains,
    )
end

function nuts_sampling_multi_threads(
    model,
    n_samples::Int,
    n_chains::Int;
    target_acceptance::Float64=0.65,
    progress=false,
)
    sampler = NUTS(target_acceptance)
    sample(model, sampler, MCMCThreads(), n_samples, n_chains, progress=progress)
end
