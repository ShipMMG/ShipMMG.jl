
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

basic_params = Mmg3DofBasicParams()
basic_params.L_pp = L_pp  # 船長Lpp[m]
basic_params.B = B  # 船幅[m]
basic_params.d = d  # 喫水[m]
basic_params.x_G = x_G  # 重心位置[]
basic_params.D_p = D_p  # プロペラ直径[m]
basic_params.m = ρ * nabla  # 質量(無次元化)[kg]
basic_params.I_zG = ρ * nabla * ((0.25 * L_pp)^2)  # 慣性モーメント[-]
basic_params.A_R = A_R  # 船の断面に対する舵面積比[-]
basic_params.η = D_p / H_R  # プロペラ直径に対する舵高さ(Dp/H)
basic_params.m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash  # 付加質量x(無次元)
basic_params.m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash  # 付加質量y(無次元)
basic_params.J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash  # 付加質量Izz(無次元)
basic_params.f_α = f_α # 直圧力勾配係数
basic_params.ϵ = ϵ  # プロペラ・舵位置伴流係数比
basic_params.t_R = t_R  # 操縦抵抗減少率
basic_params.a_H = a_H  # 舵力増加係数
basic_params.x_H = x_H_dash * L_pp  # 舵力増分作用位置
basic_params.γ_R_minus = γ_R_minus  # 整流係数
basic_params.γ_R_plus = γ_R_plus  # 整流係数
basic_params.l_R = l_r_dash  # 船長に対する舵位置
basic_params.κ = κ  # 修正係数
basic_params.t_P = t_P  # 推力減少率
basic_params.w_P0 = w_P0  # 有効伴流率
basic_params.x_P = x_P_dash  # 船長に対するプロペラ位置

maneuvering_params = Mmg3DofManeuveringParams()
maneuvering_params.k_0 = 0.2931
maneuvering_params.k_1 = -0.2753
maneuvering_params.k_2 = -0.1385
maneuvering_params.R_0_dash = 0.022
maneuvering_params.X_vv_dash = -0.040
maneuvering_params.X_vr_dash = 0.002
maneuvering_params.X_rr_dash = 0.011
maneuvering_params.X_vvvv_dash = 0.771
maneuvering_params.Y_v_dash = -0.315
maneuvering_params.Y_r_dash = 0.083
maneuvering_params.Y_vvv_dash = -1.607
maneuvering_params.Y_vvr_dash = 0.379
maneuvering_params.Y_vrr_dash = -0.391
maneuvering_params.Y_rrr_dash = 0.008
maneuvering_params.N_v_dash = -0.137
maneuvering_params.N_r_dash = -0.049
maneuvering_params.N_vvv_dash = -0.030
maneuvering_params.N_vvr_dash = -0.294
maneuvering_params.N_vrr_dash = 0.055
maneuvering_params.N_rrr_dash = -0.013
# ------

@testset "mmg.jl KVLCC2_L7 turning" begin

    duration = 200  # [s]
    max_δ_rad = 35 * pi / 180.0  # [rad]
    n_const = 17.95  # [rpm]

    sampling = duration * 10
    time_list = range(0.00, stop = duration, length = sampling)
    δ_rad_list = max_δ_rad .* ones(Float64, sampling)
    npm_list = n_const .* ones(Float64, sampling)
    mmg_results = mmg_3dof_simulate(
        time_list,
        npm_list,
        δ_rad_list,
        basic_params,
        maneuvering_params,
        u0 = 2.29 * 0.512,
        v0 = 0.0,
        r0 = 0.0,
    )
    time, u, v, r, δ, npm = mmg_results
    x, y, ψ = calc_position(time, u, v, r)
end

@testset "mmg_zigzag_test" begin
    target_δ_rad = 20.0 * π / 180.0
    target_ψ_rad_deviation = 20.0 * π / 180.0
    start_time_second = 0.00
    time_second_interval = 0.01
    end_time_second = 80.00
    time_list = start_time_second:time_second_interval:end_time_second
    n_const = 17.95  # [rpm]
    npm_list = n_const * ones(Float64, length(time_list))
    δ_list, u_list, v_list, r_list, ψ_list = mmg_3dof_zigzag_test(
        basic_params,
        maneuvering_params,
        npm_list,
        target_δ_rad,
        target_ψ_rad_deviation,
        time_second_interval,
        end_time_second,
    )
end

