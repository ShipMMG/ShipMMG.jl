"""
    mmg_3dof_model!(dX, X, p, t)

MMG 3DOF model on DifferentialEquations.ODEProblem. Update `dX`.

# Arguments
- `dX`: [du, dv, dr, dδ, dnpm]
- `X`: the initial state values. [`u`, `v`, `r`, `δ`, `npm`].
- `p`: ρ and the basic & maneuvering parameters and δ & npm spline info.
    - ρ
    - L_pp
    - B
    - d
    - x_G
    - D_p
    - m
    - I_zG
    - A_R
    - η
    - m_x
    - m_y
    - J_z
    - f_α
    - ϵ
    - t_R
    - a_H
    - x_H
    - γ_R_minus
    - γ_R_plus
    - l_R
    - κ
    - t_P
    - w_P0
    - x_P
    - k_0
    - k_1
    - k_2
    - R_0_dash
    - X_vv_dash
    - X_vr_dash
    - X_rr_dash
    - X_vvvv_dash
    - Y_v_dash
    - Y_r_dash
    - Y_vvv_dash
    - Y_vvr_dash
    - Y_vrr_dash
    - Y_rrr_dash
    - N_v_dash
    - N_r_dash
    - N_vvv_dash
    - N_vvr_dash
    - N_vrr_dash
    - N_rrr_dash
    - spl_δ
    - spl_npm
- `t`: the time.
"""
function mmg_3dof_model!(dX, X, p, t)
    u, v, r, δ, npm = X
    ρ,
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
    spl_δ,
    spl_npm = p

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

"""
Basic parameters of target ship for MMG 3DOF simulation.

# Arguments
- `L_pp::Float64`: L_pp
- `B::Float64`: 
- `d::Float64`:
- `x_G::Float64`: 
- `D_p::Float64`: 
- `m::Float64`: 
- `I_zG::Float64`: 
- `A_R::Float64`: 
- `η::Float64`: 
- `m_x::Float64`: 
- `m_y::Float64`: 
- `J_z::Float64`: 
- `f_α::Float64`: 
- `ϵ::Float64`: 
- `t_R::Float64`: 
- `a_H::Float64`: 
- `x_H::Float64`: 
- `γ_R_minus::Float64`: 
- `γ_R_plus::Float64`: 
- `l_R::Float64`: 
- `κ::Float64`: 
- `t_P::Float64`: 
- `w_P0::Float64`: 
- `x_P::Float64`: 
"""
mutable struct Mmg3DofBasicParams
    L_pp::Float64
    B::Float64
    d::Float64
    x_G::Float64
    D_p::Float64
    m::Float64
    I_zG::Float64
    A_R::Float64
    η::Float64
    m_x::Float64
    m_y::Float64
    J_z::Float64
    f_α::Float64
    ϵ::Float64
    t_R::Float64
    a_H::Float64
    x_H::Float64
    γ_R_minus::Float64
    γ_R_plus::Float64
    l_R::Float64
    κ::Float64
    t_P::Float64
    w_P0::Float64
    x_P::Float64
    Mmg3DofBasicParams() = new()
end

"""
Maneuvering parameters of target ship for MMG 3DOF simulation.

# Arguments
- `k_0::Float64`
- `k_1::Float64`
- `k_2::Float64`
- `R_0_dash::Float64`
- `X_vv_dash::Float64`
- `X_vr_dash::Float64`
- `X_rr_dash::Float64`
- `X_vvvv_dash::Float64`
- `Y_v_dash::Float64`
- `Y_r_dash::Float64`
- `Y_vvv_dash::Float64`
- `Y_vvr_dash::Float64`
- `Y_vrr_dash::Float64`
- `Y_rrr_dash::Float64`
- `N_v_dash::Float64`
- `N_r_dash::Float64`
- `N_vvv_dash::Float64`
- `N_vvr_dash::Float64`
- `N_vrr_dash::Float64`
- `N_rrr_dash::Float64`
"""
mutable struct Mmg3DofManeuveringParams
    k_0::Float64
    k_1::Float64
    k_2::Float64
    R_0_dash::Float64
    X_vv_dash::Float64
    X_vr_dash::Float64
    X_rr_dash::Float64
    X_vvvv_dash::Float64
    Y_v_dash::Float64
    Y_r_dash::Float64
    Y_vvv_dash::Float64
    Y_vvr_dash::Float64
    Y_vrr_dash::Float64
    Y_rrr_dash::Float64
    N_v_dash::Float64
    N_r_dash::Float64
    N_vvv_dash::Float64
    N_vvr_dash::Float64
    N_vrr_dash::Float64
    N_rrr_dash::Float64
    Mmg3DofManeuveringParams() = new()
end

"""
    mmg_3dof_simulate(time_list, npm_list, δ_list, basic_params, maneuvering_params, [, u0, v0, r0, ρ, algorithm, reltol, abstol]) -> time, u, v, r, δ, npm

Returns the MMG 3DOF simulation results including the lists of time, u, v, r, δ, npm.
This function has the same logic of `ShipMMG.simulate()`.

# Arguments
- `basic_params::Mmg3DofBasicParams`: the basic parameters of target ship.
- `maneuvering_params::Mmg3DofManeuveringParams`: the maneuvering parameters of target ship.
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `npm_list`: the list of propeller rpm.
- `u0::Float64=0.0`: the initial x (surge) velocity.
- `v0::Float64=0.0`: the initial y (sway) velocity.
- `r0::Float64=0.0`: the initial rate of turn [rad/s].
- `ρ::Float64=1.025`: the seawater density [g/cm^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KVLCC2_L7 turning test.

```julia-rep1
julia> ρ = 1.025;
julia> L_pp = 7.00;
julia> B = 1.27;
julia> d = 0.46;
julia> nabla = 3.27;
julia> x_G = 0.25;
julia> # C_b = 0.810;
julia> D_p = 0.216;
julia> H_R = 0.345;
julia> A_R = 0.0539;
julia> t_P = 0.220;
julia> w_P0 = 0.40;
julia> m_x_dash = 0.022;
julia> m_y_dash = 0.223;
julia> J_z_dash = 0.011;
julia> t_R = 0.387;
julia> a_H = 0.312;
julia> x_H_dash = -0.464;
julia> γ_R_minus = 0.395;
julia> γ_R_plus = 0.640;
julia> l_r_dash = -0.710;
julia> x_P_dash = -0.690;
julia> ϵ = 1.09;
julia> κ = 0.50;
julia> f_α = 2.747;
julia> basic_params = Mmg3DofBasicParams();
julia> basic_params.L_pp = L_pp;  # 船長Lpp[m]
julia> basic_params.B = B;  # 船幅[m]
julia> basic_params.d = d;  # 喫水[m]
julia> basic_params.x_G = x_G;  # 重心位置[]
julia> basic_params.D_p = D_p;  # プロペラ直径[m]
julia> basic_params.m = ρ * nabla;  # 質量(無次元化)[kg]
julia> basic_params.I_zG = ρ * nabla * ((0.25 * L_pp)^2);  # 慣性モーメント[-]
julia> basic_params.A_R = A_R;  # 船の断面に対する舵面積比[-]
julia> basic_params.η = D_p / H_R;  # プロペラ直径に対する舵高さ(Dp/H)
julia> basic_params.m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash;  # 付加質量x(無次元)
julia> basic_params.m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash;  # 付加質量y(無次元)
julia> basic_params.J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash;  # 付加質量Izz(無次元)
julia> basic_params.f_α = f_α; # 直圧力勾配係数
julia> basic_params.ϵ = ϵ;  # プロペラ・舵位置伴流係数比
julia> basic_params.t_R = t_R;  # 操縦抵抗減少率
julia> basic_params.a_H = a_H;  # 舵力増加係数
julia> basic_params.x_H = x_H_dash * L_pp;  # 舵力増分作用位置
julia> basic_params.γ_R_minus = γ_R_minus;  # 整流係数
julia> basic_params.γ_R_plus = γ_R_plus;  # 整流係数
julia> basic_params.l_R = l_r_dash;  # 船長に対する舵位置
julia> basic_params.κ = κ;  # 修正係数
julia> basic_params.t_P = t_P;  # 推力減少率
julia> basic_params.w_P0 = w_P0;  # 有効伴流率
julia> basic_params.x_P = x_P_dash;  # 船長に対するプロペラ位置
julia> maneuvering_params = Mmg3DofManeuveringParams();
julia> maneuvering_params.k_0 = 0.2931;
julia> maneuvering_params.k_1 = -0.2753;
julia> maneuvering_params.k_2 = -0.1385;
julia> maneuvering_params.R_0_dash = 0.022;
julia> maneuvering_params.X_vv_dash = -0.040;
julia> maneuvering_params.X_vr_dash = 0.002;
julia> maneuvering_params.X_rr_dash = 0.011;
julia> maneuvering_params.X_vvvv_dash = 0.771;
julia> maneuvering_params.Y_v_dash = -0.315;
julia> maneuvering_params.Y_r_dash = 0.083;
julia> maneuvering_params.Y_vvv_dash = -1.607;
julia> maneuvering_params.Y_vvr_dash = 0.379;
julia> maneuvering_params.Y_vrr_dash = -0.391;
julia> maneuvering_params.Y_rrr_dash = 0.008;
julia> maneuvering_params.N_v_dash = -0.137;
julia> maneuvering_params.N_r_dash = -0.049;
julia> maneuvering_params.N_vvv_dash = -0.030;
julia> maneuvering_params.N_vvr_dash = -0.294;
julia> maneuvering_params.N_vrr_dash = 0.055;
julia> maneuvering_params.N_rrr_dash = -0.013;
julia> duration = 200; # [s]
julia> max_δ_rad = 35 * pi / 180.0;  # [rad]
julia> n_const = 17.95;  # [rpm]
julia> sampling = duration * 10;
julia> time_list = range(0.00, stop = duration, length = sampling);
julia> δ_rad_list = max_δ_rad .* ones(Float64, sampling);
julia> npm_list = n_const .* ones(Float64, sampling);
julia> mmg_results = mmg_3dof_simulate(
    basic_params,
    maneuvering_params,
    time_list,
    δ_rad_list,
    npm_list,
    u0 = 2.29 * 0.512,
    v0 = 0.0,
    r0 = 0.0,
);
```
"""
function mmg_3dof_simulate(
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams,
    time_list,
    δ_list,
    npm_list;
    u0::Float64 = 0.0,
    v0::Float64 = 0.0,
    r0::Float64 = 0.0,
    ρ::Float64 = 1.025,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)
    simulate(
        basic_params.L_pp,
        basic_params.B,
        basic_params.d,
        basic_params.x_G,
        basic_params.D_p,
        basic_params.m,
        basic_params.I_zG,
        basic_params.A_R,
        basic_params.η,
        basic_params.m_x,
        basic_params.m_y,
        basic_params.J_z,
        basic_params.f_α,
        basic_params.ϵ,
        basic_params.t_R,
        basic_params.a_H,
        basic_params.x_H,
        basic_params.γ_R_minus,
        basic_params.γ_R_plus,
        basic_params.l_R,
        basic_params.κ,
        basic_params.t_P,
        basic_params.w_P0,
        basic_params.x_P,
        maneuvering_params.k_0,
        maneuvering_params.k_1,
        maneuvering_params.k_2,
        maneuvering_params.R_0_dash,
        maneuvering_params.X_vv_dash,
        maneuvering_params.X_vr_dash,
        maneuvering_params.X_rr_dash,
        maneuvering_params.X_vvvv_dash,
        maneuvering_params.Y_v_dash,
        maneuvering_params.Y_r_dash,
        maneuvering_params.Y_vvv_dash,
        maneuvering_params.Y_vvr_dash,
        maneuvering_params.Y_vrr_dash,
        maneuvering_params.Y_rrr_dash,
        maneuvering_params.N_v_dash,
        maneuvering_params.N_r_dash,
        maneuvering_params.N_vvv_dash,
        maneuvering_params.N_vvr_dash,
        maneuvering_params.N_vrr_dash,
        maneuvering_params.N_rrr_dash,
        time_list,
        δ_list,
        npm_list,
        u0 = u0,
        v0 = v0,
        r0 = r0,
        ρ = ρ,
        algorithm = algorithm,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    mmg_3dof_simulate(time_list, npm_list, δ_list, L_pp, B, d, x_G, D_p, m, I_zG, A_R, η, m_x, m_y, J_z, f_α, ϵ, t_R, a_H, x_H, γ_R_minus, γ_R_plus, l_R, κ, t_P, w_P0, x_P, k_0, k_1, k_2, R_0_dash, X_vv_dash, X_vr_dash, X_rr_dash, X_vvvv_dash, Y_v_dash, Y_r_dash, Y_vvv_dash, Y_vvr_dash, Y_vrr_dash, Y_rrr_dash, N_v_dash, N_r_dash, N_vvv_dash, N_vvr_dash, N_vrr_dash, N_rrr_dash, [, u0, v0, r0, ρ, algorithm, reltol, abstol]) -> time, u, v, r, δ, npm

Returns the MMG 3DOF simulation results including the lists of time, u, v, r, δ, npm.
This function has the same logic of `ShipMMG.mmg_3dof_simulate()`.

# Arguments
- `L_pp::Float64`: L_pp
- `B::Float64`: 
- `d::Float64`:
- `x_G::Float64`: 
- `D_p::Float64`: 
- `m::Float64`: 
- `I_zG::Float64`: 
- `A_R::Float64`: 
- `η::Float64`: 
- `m_x::Float64`: 
- `m_y::Float64`: 
- `J_z::Float64`: 
- `f_α::Float64`: 
- `ϵ::Float64`: 
- `t_R::Float64`: 
- `a_H::Float64`: 
- `x_H::Float64`: 
- `γ_R_minus::Float64`: 
- `γ_R_plus::Float64`: 
- `l_R::Float64`: 
- `κ::Float64`: 
- `t_P::Float64`: 
- `w_P0::Float64`: 
- `x_P::Float64`: 
- `k_0::Float64`
- `k_1::Float64`
- `k_2::Float64`
- `R_0_dash::Float64`
- `X_vv_dash::Float64`
- `X_vr_dash::Float64`
- `X_rr_dash::Float64`
- `X_vvvv_dash::Float64`
- `Y_v_dash::Float64`
- `Y_r_dash::Float64`
- `Y_vvv_dash::Float64`
- `Y_vvr_dash::Float64`
- `Y_vrr_dash::Float64`
- `Y_rrr_dash::Float64`
- `N_v_dash::Float64`
- `N_r_dash::Float64`
- `N_vvv_dash::Float64`
- `N_vvr_dash::Float64`
- `N_vrr_dash::Float64`
- `N_rrr_dash::Float64`
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `npm_list`: the list of propeller rpm.
- `u0::Float64=0.0`: the initial x (surge) velocity.
- `v0::Float64=0.0`: the initial y (sway) velocity.
- `r0::Float64=0.0`: the initial rate of turn [rad/s].
- `ρ::Float64=1.025`: the seawater density [g/cm^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
"""
function simulate(
    L_pp::Float64,
    B::Float64,
    d::Float64,
    x_G::Float64,
    D_p::Float64,
    m::Float64,
    I_zG::Float64,
    A_R::Float64,
    η::Float64,
    m_x::Float64,
    m_y::Float64,
    J_z::Float64,
    f_α::Float64,
    ϵ::Float64,
    t_R::Float64,
    a_H::Float64,
    x_H::Float64,
    γ_R_minus::Float64,
    γ_R_plus::Float64,
    l_R::Float64,
    κ::Float64,
    t_P::Float64,
    w_P0::Float64,
    x_P::Float64,
    k_0::Float64,
    k_1::Float64,
    k_2::Float64,
    R_0_dash::Float64,
    X_vv_dash::Float64,
    X_vr_dash::Float64,
    X_rr_dash::Float64,
    X_vvvv_dash::Float64,
    Y_v_dash::Float64,
    Y_r_dash::Float64,
    Y_vvv_dash::Float64,
    Y_vvr_dash::Float64,
    Y_vrr_dash::Float64,
    Y_rrr_dash::Float64,
    N_v_dash::Float64,
    N_r_dash::Float64,
    N_vvv_dash::Float64,
    N_vvr_dash::Float64,
    N_vrr_dash::Float64,
    N_rrr_dash::Float64,
    time_list,
    δ_list,
    npm_list;
    u0::Float64 = 0.0,
    v0::Float64 = 0.0,
    r0::Float64 = 0.0,
    ρ::Float64 = 1.025,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)
    spl_δ = Spline1D(time_list, δ_list)
    spl_npm = Spline1D(time_list, npm_list)

    X0 = [u0; v0; r0; δ_list[1]; npm_list[1]]
    p = [
        ρ,
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
        spl_δ,
        spl_npm,
    ]
    prob = ODEProblem(mmg_3dof_model!, X0, (time_list[1], time_list[end]), p)
    sol = solve(
        prob,
        algorithm,
        reltol = reltol,
        abstol = abstol,
        saveat = time_list[2] - time_list[1],
    )
    results = hcat(sol.u...)
    time = sol.t
    u = results[1, :]
    v = results[2, :]
    r = results[3, :]
    δ = results[4, :]
    npm = results[5, :]
    time, u, v, r, δ, npm
end

"""
    mmg_3dof_zigzag_test(basic_params, maneuvering_params, time_list, npm_list, target_δ_rad, target_ψ_rad_deviation, [, u0, v0, r0, ψ0, δ0, δ_rad_rate, algorithm, reltol, abstol]) -> final_time_list, final_u_list, final_v_list, final_r_list, final_ψ_list, final_δ_list

Returns the MMG 3DOF zigzag simulation results.

# Arguments
- `basic_params::Mmg3DofBasicParams`: the basic parameters of target ship.
- `maneuvering_params::Mmg3DofManeuveringParams`: the maneuvering parameters of target ship.
- `time_list`: the list of simulatino time.
- `npm_list`: the list of propeller rpm.
- `target_δ_rad::Float64`: target rudder angle of zigzag test.
- `target_ψ_rad_deviation::Float64`: target azimuth deviation of zigzag test.
- `u0::Float64=0.0`: the initial x (surge) velocity.
- `v0::Float64=0.0`: the initial y (sway) velocity.
- `r0::Float64=0.0`: the initial rate of turn [rad/s].
- `δ0::Float64=0.0`: the initial rudder angle.
- `δ_rad_rate::Float64=10.0*π/180`: the change rate of rudder angle [rad/s]. 
- `ρ::Float64=1.025`: the seawater density [g/cm^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KVLCC2_L7 zigzag test.

```julia-rep1
julia> ρ = 1.025;
julia> L_pp = 7.00;
julia> B = 1.27;
julia> d = 0.46;
julia> nabla = 3.27;
julia> x_G = 0.25;
julia> # C_b = 0.810;
julia> D_p = 0.216;
julia> H_R = 0.345;
julia> A_R = 0.0539;
julia> t_P = 0.220;
julia> w_P0 = 0.40;
julia> m_x_dash = 0.022;
julia> m_y_dash = 0.223;
julia> J_z_dash = 0.011;
julia> t_R = 0.387;
julia> a_H = 0.312;
julia> x_H_dash = -0.464;
julia> γ_R_minus = 0.395;
julia> γ_R_plus = 0.640;
julia> l_r_dash = -0.710;
julia> x_P_dash = -0.690;
julia> ϵ = 1.09;
julia> κ = 0.50;
julia> f_α = 2.747;
julia> basic_params = Mmg3DofBasicParams();
julia> basic_params.L_pp = L_pp;  # 船長Lpp[m]
julia> basic_params.B = B;  # 船幅[m]
julia> basic_params.d = d;  # 喫水[m]
julia> basic_params.x_G = x_G;  # 重心位置[]
julia> basic_params.D_p = D_p;  # プロペラ直径[m]
julia> basic_params.m = ρ * nabla;  # 質量(無次元化)[kg]
julia> basic_params.I_zG = ρ * nabla * ((0.25 * L_pp)^2);  # 慣性モーメント[-]
julia> basic_params.A_R = A_R;  # 船の断面に対する舵面積比[-]
julia> basic_params.η = D_p / H_R;  # プロペラ直径に対する舵高さ(Dp/H)
julia> basic_params.m_x = (0.5 * ρ * (L_pp^2) * d) * m_x_dash;  # 付加質量x(無次元)
julia> basic_params.m_y = (0.5 * ρ * (L_pp^2) * d) * m_y_dash;  # 付加質量y(無次元)
julia> basic_params.J_z = (0.5 * ρ * (L_pp^4) * d) * J_z_dash;  # 付加質量Izz(無次元)
julia> basic_params.f_α = f_α; # 直圧力勾配係数
julia> basic_params.ϵ = ϵ;  # プロペラ・舵位置伴流係数比
julia> basic_params.t_R = t_R;  # 操縦抵抗減少率
julia> basic_params.a_H = a_H;  # 舵力増加係数
julia> basic_params.x_H = x_H_dash * L_pp;  # 舵力増分作用位置
julia> basic_params.γ_R_minus = γ_R_minus;  # 整流係数
julia> basic_params.γ_R_plus = γ_R_plus;  # 整流係数
julia> basic_params.l_R = l_r_dash;  # 船長に対する舵位置
julia> basic_params.κ = κ;  # 修正係数
julia> basic_params.t_P = t_P;  # 推力減少率
julia> basic_params.w_P0 = w_P0;  # 有効伴流率
julia> basic_params.x_P = x_P_dash;  # 船長に対するプロペラ位置
julia> maneuvering_params = Mmg3DofManeuveringParams();
julia> maneuvering_params.k_0 = 0.2931;
julia> maneuvering_params.k_1 = -0.2753;
julia> maneuvering_params.k_2 = -0.1385;
julia> maneuvering_params.R_0_dash = 0.022;
julia> maneuvering_params.X_vv_dash = -0.040;
julia> maneuvering_params.X_vr_dash = 0.002;
julia> maneuvering_params.X_rr_dash = 0.011;
julia> maneuvering_params.X_vvvv_dash = 0.771;
julia> maneuvering_params.Y_v_dash = -0.315;
julia> maneuvering_params.Y_r_dash = 0.083;
julia> maneuvering_params.Y_vvv_dash = -1.607;
julia> maneuvering_params.Y_vvr_dash = 0.379;
julia> maneuvering_params.Y_vrr_dash = -0.391;
julia> maneuvering_params.Y_rrr_dash = 0.008;
julia> maneuvering_params.N_v_dash = -0.137;
julia> maneuvering_params.N_r_dash = -0.049;
julia> maneuvering_params.N_vvv_dash = -0.030;
julia> maneuvering_params.N_vvr_dash = -0.294;
julia> maneuvering_params.N_vrr_dash = 0.055;
julia> maneuvering_params.N_rrr_dash = -0.013;
julia> target_δ_rad = 20.0 * π / 180.0
julia> target_ψ_rad_deviation = 20.0 * π / 180.0
julia> start_time_second = 0.00
julia> time_second_interval = 0.01
julia> end_time_second = 80.00
julia> time_list = start_time_second:time_second_interval:end_time_second
julia> n_const = 17.95  # [rpm]
julia> npm_list = n_const * ones(Float64, length(time_list))
julia> time_list, δ_list, u_list, v_list, r_list, ψ_list = mmg_3dof_zigzag_test(
    basic_params,
    maneuvering_params,
    time_list
    npm_list,
    target_δ_rad,
    target_ψ_rad_deviation,
);
```
"""
function mmg_3dof_zigzag_test(
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams,
    time_list,
    npm_list,
    target_δ_rad::Float64,
    target_ψ_rad_deviation::Float64;
    u0::Float64 = 0.0,
    v0::Float64 = 0.0,
    r0::Float64 = 0.0,
    ψ0::Float64 = 0.0,
    δ0::Float64 = 0.0,
    δ_rad_rate::Float64 = 10.0 * π / 180,
    ρ::Float64 = 1.025,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)
    target_ψ_rad_deviation = abs(target_ψ_rad_deviation)

    # time_list = start_time_second:time_second_interval:end_time_second
    final_time_list = zeros(length(time_list))
    final_δ_list = zeros(length(time_list))
    final_u_list = zeros(length(time_list))
    final_v_list = zeros(length(time_list))
    final_r_list = zeros(length(time_list))
    final_ψ_list = zeros(length(time_list))

    next_stage_index = 1
    target_δ_rad = -target_δ_rad  # for changing in while loop
    ψ = ψ0
    while next_stage_index < length(time_list)
        target_δ_rad = -target_δ_rad
        start_index = next_stage_index

        # Make delta list
        δ_list = zeros(length(time_list) - start_index + 1)
        if start_index == 1
            δ_list[1] = δ0
            u0 = u0
            v0 = v0
            r0 = r0
        else
            δ_list[1] = final_δ_list[start_index-1]
            u0 = final_u_list[start_index-1]
            v0 = final_v_list[start_index-1]
            r0 = final_r_list[start_index-1]
        end

        for i = (start_index+1):length(time_list)
            Δt = time_list[i] - time_list[i-1]
            if target_δ_rad > 0
                δ = δ_list[i-start_index] + δ_rad_rate * Δt
                if δ >= target_δ_rad
                    δ = target_δ_rad
                end
                δ_list[i-start_index+1] = δ
            elseif target_δ_rad <= 0
                δ = δ_list[i-start_index] - δ_rad_rate * Δt
                if δ <= target_δ_rad
                    δ = target_δ_rad
                end
                δ_list[i-start_index+1] = δ
            end
        end

        time, u, v, r, δ, npm = mmg_3dof_simulate(
            basic_params,
            maneuvering_params,
            time_list[start_index:end],
            δ_list,
            npm_list[start_index:end],
            u0 = u0,
            v0 = v0,
            r0 = r0,
            ρ = ρ,
            algorithm = algorithm,
            reltol = reltol,
            abstol = abstol,
        )
        x, y, ψ_list = calc_position(time, u, v, r, x0 = 0.0, y0 = 0.0, ψ0 = ψ)
        # get finish index
        target_ψ_rad = ψ0 + target_ψ_rad_deviation
        if target_δ_rad < 0
            target_ψ_rad = ψ0 - target_ψ_rad_deviation
        end
        over_index = findfirst(e -> e >= target_ψ_rad, ψ_list)
        if target_δ_rad < 0
            over_index = findfirst(e -> e <= target_ψ_rad, ψ_list)
        end
        next_stage_index = length(time_list)
        if isnothing(over_index)
            final_time_list[start_index:next_stage_index] = time
            final_δ_list[start_index:next_stage_index] = δ_list
            final_u_list[start_index:next_stage_index] = u
            final_v_list[start_index:next_stage_index] = v
            final_r_list[start_index:next_stage_index] = r
            final_ψ_list[start_index:next_stage_index] = ψ_list
        else
            ψ = ψ_list[over_index]
            next_stage_index = over_index + start_index - 1
            final_time_list[start_index:next_stage_index] = time[begin:over_index]
            final_δ_list[start_index:next_stage_index] = δ_list[begin:over_index]
            final_u_list[start_index:next_stage_index] = u[begin:over_index]
            final_v_list[start_index:next_stage_index] = v[begin:over_index]
            final_r_list[start_index:next_stage_index] = r[begin:over_index]
            final_ψ_list[start_index:next_stage_index] = ψ_list[begin:over_index]
        end
    end
    final_time_list, final_u_list, final_v_list, final_r_list, final_ψ_list, final_δ_list
end
