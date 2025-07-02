"""
    mmg_3dof_model!(dX, X, p, t)

MMG 3DOF model on DifferentialEquations.ODEProblem. Update `dX`.

# Arguments
- `dX`: [du, dv, dr, dx, dy, dΨ, dδ, dn_p]
- `X`: the initial state values. [`u`, `v`, `r`, `x`, `y`, `Ψ`, `δ`, `n_p`].
- `p`: ρ and the basic & maneuvering parameters and δ & n_p spline info.
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
    - x_R
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
    - spl_n_p
- `t`: the time.
"""
function mmg_3dof_model!(dX, X, p, t)
    u, v, r, x, y, Ψ, δ, n_p = X
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
    spl_n_p = p

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
    if n_p == 0.0
        J = 0.0
    else
        J = (1 - w_P) * u / (n_p * D_p)
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
        u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * n_p^2 * D_p^4 / pi)^2)
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
    X_P = (1.0 - t_P) * ρ * K_T * n_p^2 * D_p^4
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
    N_R = -(x_R + a_H * x_H) * F_N * cos(δ)
    dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
    dX[2] =
        dv =
            (
                (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
            ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
    dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
    dX[4] = dx = u * cos(Ψ) - v * sin(Ψ)
    dX[5] = dy = u * sin(Ψ) + v * cos(Ψ)
    dX[6] = dΨ = r
    dX[7] = dδ = derivative(spl_δ, t)
    dX[8] = dn_p = derivative(spl_n_p, t)
end

"""
Basic parameters of target ship for MMG 3DOF simulation.

# Arguments
- `L_pp`: L_pp
- `B`: 
- `d`:
- `x_G`: 
- `D_p`: 
- `m`: 
- `I_zG`: 
- `A_R`: 
- `η`: 
- `m_x`: 
- `m_y`: 
- `J_z`: 
- `f_α`: 
- `ϵ`: 
- `t_R`: 
- `x_R`: 
- `a_H`: 
- `x_H`: 
- `γ_R_minus`: 
- `γ_R_plus`: 
- `l_R`: 
- `κ`: 
- `t_P`: 
- `w_P0`: 
- `x_P`: 
"""
@with_kw mutable struct Mmg3DofBasicParams{T}
    L_pp::T
    B::T
    d::T
    x_G::T
    D_p::T
    m::T
    I_zG::T
    A_R::T
    η::T
    m_x::T
    m_y::T
    J_z::T
    f_α::T
    ϵ::T
    t_R::T
    x_R::T
    a_H::T
    x_H::T
    γ_R_minus::T
    γ_R_plus::T
    l_R::T
    κ::T
    t_P::T
    w_P0::T
    x_P::T
end

"""
Maneuvering parameters of target ship for MMG 3DOF simulation.

# Arguments
- `k_0::T`
- `k_1::T`
- `k_2::T`
- `R_0_dash::T`
- `X_vv_dash::T`
- `X_vr_dash::T`
- `X_rr_dash::T`
- `X_vvvv_dash::T`
- `Y_v_dash::T`
- `Y_r_dash::T`
- `Y_vvv_dash::T`
- `Y_vvr_dash::T`
- `Y_vrr_dash::T`
- `Y_rrr_dash::T`
- `N_v_dash::T`
- `N_r_dash::T`
- `N_vvv_dash::T`
- `N_vvr_dash::T`
- `N_vrr_dash::T`
- `N_rrr_dash::T`
"""
@with_kw mutable struct Mmg3DofManeuveringParams{T}
    k_0::T
    k_1::T
    k_2::T
    R_0_dash::T
    X_vv_dash::T
    X_vr_dash::T
    X_rr_dash::T
    X_vvvv_dash::T
    Y_v_dash::T
    Y_r_dash::T
    Y_vvv_dash::T
    Y_vvr_dash::T
    Y_vrr_dash::T
    Y_rrr_dash::T
    N_v_dash::T
    N_r_dash::T
    N_vvv_dash::T
    N_vvr_dash::T
    N_vrr_dash::T
    N_rrr_dash::T
end

"""
    mmg_3dof_simulate(time_list, n_p_list, δ_list, basic_params, maneuvering_params, [, u0, v0, r0, x0, y0, Ψ0, ρ, algorithm, reltol, abstol]) -> u, v, r, x, y, Ψ, δ, n_p

Returns the MMG 3DOF simulation results including the lists of time, u, v, r, x, y, Ψ, δ, n_p.
This function has the same logic of `ShipMMG.simulate()`.

# Arguments
- `basic_params::Mmg3DofBasicParams`: the basic parameters of target ship.
- `maneuvering_params::Mmg3DofManeuveringParams`: the maneuvering parameters of target ship.
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `n_p_list`: the list of propeller rps.
- `u0=0.0`: the initial x (surge) velocity.
- `v0=0.0`: the initial y (sway) velocity.
- `r0=0.0`: the initial rate of turn [rad/s].
- `x0=0.0`: the initial x (surge) position.
- `y0=0.0`: the initial y (sway) position.
- `Ψ0=0.0`: the initial Ψ (yaw) azimuth [rad].
- `ρ=1025.0`: the seawater density [kg/m^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KVLCC2_L7 turning test.

```julia-rep1
julia> basic_params, maneuvering_params = get_KVLCC2_L7_params();
julia> duration = 200; # [s]
julia> max_δ_rad = 35 * pi / 180.0;  # [rad]
julia> n_const = 17.95;  # [rps]
julia> sampling = duration * 10;
julia> time_list = range(0.00, stop = duration, length = sampling);
julia> δ_rad_list = max_δ_rad .* ones(Float64, sampling);
julia> n_p_list = n_const .* ones(Float64, sampling);
julia> mmg_results = mmg_3dof_simulate(
    basic_params,
    maneuvering_params,
    time_list,
    δ_rad_list,
    n_p_list,
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
    n_p_list;
    u0=0.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    ρ=1025.0,
    algorithm=Tsit5(),
    reltol=1e-8,
    abstol=1e-8
)
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
    x_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    @unpack k_0,
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
    N_rrr_dash = maneuvering_params
    simulate(
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
        time_list,
        δ_list,
        n_p_list,
        u0=u0,
        v0=v0,
        r0=r0,
        x0=x0,
        y0=y0,
        Ψ0=Ψ0,
        ρ=ρ,
        algorithm=algorithm,
        reltol=reltol,
        abstol=abstol,
    )
end

"""
    simulate(time_list, n_p_list, δ_list, L_pp, B, d, x_G, D_p, m, I_zG, A_R, η, m_x, m_y, J_z, f_α, ϵ, t_R, x_R, a_H, x_H, γ_R_minus, γ_R_plus, l_R, κ, t_P, w_P0, x_P, k_0, k_1, k_2, R_0_dash, X_vv_dash, X_vr_dash, X_rr_dash, X_vvvv_dash, Y_v_dash, Y_r_dash, Y_vvv_dash, Y_vvr_dash, Y_vrr_dash, Y_rrr_dash, N_v_dash, N_r_dash, N_vvv_dash, N_vvr_dash, N_vrr_dash, N_rrr_dash, [, u0, v0, r0, ρ, algorithm, reltol, abstol]) -> u, v, r, x, y, Ψ, δ, n_p

Returns the MMG 3DOF simulation results including the lists of time, u, v, r, x, y, Ψ, δ, n_p.
This function has the same logic of `ShipMMG.mmg_3dof_simulate()`.

# Arguments
- `L_pp`: L_pp
- `B`: 
- `d`:
- `x_G`: 
- `D_p`: 
- `m`: 
- `I_zG`: 
- `A_R`: 
- `η`: 
- `m_x`: 
- `m_y`: 
- `J_z`: 
- `f_α`: 
- `ϵ`: 
- `t_R`: 
- `x_R`: 
- `a_H`: 
- `x_H`: 
- `γ_R_minus`: 
- `γ_R_plus`: 
- `l_R`: 
- `κ`: 
- `t_P`: 
- `w_P0`: 
- `x_P`: 
- `k_0`
- `k_1`
- `k_2`
- `R_0_dash`
- `X_vv_dash`
- `X_vr_dash`
- `X_rr_dash`
- `X_vvvv_dash`
- `Y_v_dash`
- `Y_r_dash`
- `Y_vvv_dash`
- `Y_vvr_dash`
- `Y_vrr_dash`
- `Y_rrr_dash`
- `N_v_dash`
- `N_r_dash`
- `N_vvv_dash`
- `N_vvr_dash`
- `N_vrr_dash`
- `N_rrr_dash`
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `n_p_list`: the list of propeller rps.
- `u0=0.0`: the initial x (surge) velocity.
- `v0=0.0`: the initial y (sway) velocity.
- `r0=0.0`: the initial rate of turn [rad/s].
- `x0=0.0`: the initial x (surge) position.
- `y0=0.0`: the initial y (sway) position.
- `Ψ0=0.0`: the initial Ψ (yaw) azimuth [rad].
- `ρ=1025.0`: the seawater density [kg/m^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
"""
function simulate(
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
    time_list,
    δ_list,
    n_p_list;
    u0=0.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    ρ=1025.0,
    algorithm=Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
)
    spl_δ = Spline1D(time_list, δ_list)
    spl_n_p = Spline1D(time_list, n_p_list)

    X0 = [u0; v0; r0; x0; y0; Ψ0; δ_list[1]; n_p_list[1]]
    p = (
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
        spl_n_p,
    )
    prob = ODEProblem(mmg_3dof_model!, X0, (time_list[1], time_list[end]), p)
    sol = solve(prob, algorithm, reltol=reltol, abstol=abstol)
    sol_timelist = sol(time_list)
    results = hcat(sol_timelist.u...)
    u = results[1, :]
    v = results[2, :]
    r = results[3, :]
    x = results[4, :]
    y = results[5, :]
    Ψ = results[6, :]
    δ = results[7, :]
    n_p = results[8, :]
    u, v, r, x, y, Ψ, δ, n_p
end

"""
    mmg_3dof_zigzag_test(basic_params, maneuvering_params, time_list, n_p_list, target_δ_rad, target_Ψ_rad_deviation, [, u0, v0, r0, x0, y0, Ψ0, δ0, δ_rad_rate, algorithm, reltol, abstol]) -> u, v, r, x, y, Ψ, δ

Returns the MMG 3DOF zigzag simulation results.

# Arguments
- `basic_params::Mmg3DofBasicParams`: the basic parameters of target ship.
- `maneuvering_params::Mmg3DofManeuveringParams`: the maneuvering parameters of target ship.
- `time_list`: the list of simulatino time.
- `n_p_list`: the list of propeller rps.
- `target_δ_rad`: target rudder angle of zigzag test.
- `target_Ψ_rad_deviation`: target azimuth deviation of zigzag test.
- `u0=0.0`: the initial x (surge) velocity.
- `v0=0.0`: the initial y (sway) velocity.
- `r0=0.0`: the initial rate of turn [rad/s].
- `x0=0.0`: the initial x (surge) position.
- `y0=0.0`: the initial y (sway) position.
- `Ψ0=0.0`: the initial Ψ (yaw) azimuth [rad].
- `δ0=0.0`: the initial rudder angle.
- `δ_rad_rate=10.0*π/180`: the change rate of rudder angle [rad/s]. 
- `ρ=1025.0`: the seawater density [kg/m^3].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KVLCC2_L7 zigzag test.

```julia-rep1
julia> ρ = 1025.0;
julia> basic_params, maneuvering_params = get_KVLCC2_L7_params();
julia> target_δ_rad = 20.0 * π / 180.0
julia> target_Ψ_rad_deviation = 20.0 * π / 180.0
julia> start_time_second = 0.00
julia> time_second_interval = 0.01
julia> end_time_second = 80.00
julia> time_list = start_time_second:time_second_interval:end_time_second
julia> n_const = 17.95  # [rps]
julia> n_p_list = n_const * ones(Float64, length(time_list))
julia> δ_list, u_list, v_list, r_list, Ψ_list = mmg_3dof_zigzag_test(
    basic_params,
    maneuvering_params,
    time_list
    n_p_list,
    target_δ_rad,
    target_Ψ_rad_deviation,
);
```
"""
function mmg_3dof_zigzag_test(
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams,
    time_list,
    n_p_list,
    target_δ_rad,
    target_Ψ_rad_deviation;
    u0=0.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    δ0=0.0,
    δ_rad_rate=10.0 * π / 180,
    ρ=1025.0,
    algorithm=Tsit5(),
    reltol=1e-8,
    abstol=1e-8
)
    target_Ψ_rad_deviation = abs(target_Ψ_rad_deviation)

    final_δ_list = zeros(length(time_list))
    final_u_list = zeros(length(time_list))
    final_v_list = zeros(length(time_list))
    final_r_list = zeros(length(time_list))
    final_x_list = zeros(length(time_list))
    final_y_list = zeros(length(time_list))
    final_Ψ_list = zeros(length(time_list))

    next_stage_index = 1
    target_δ_rad = -target_δ_rad  # for changing in while loop
    Ψ = Ψ0
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
            x0 = x0
            y0 = y0
        else
            δ_list[1] = final_δ_list[start_index-1]
            u0 = final_u_list[start_index-1]
            v0 = final_v_list[start_index-1]
            r0 = final_r_list[start_index-1]
            x0 = final_x_list[start_index-1]
            y0 = final_y_list[start_index-1]
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

        u, v, r, x, y, Ψ_list, δ, n_p = mmg_3dof_simulate(
            basic_params,
            maneuvering_params,
            time_list[start_index:end],
            δ_list,
            n_p_list[start_index:end],
            u0=u0,
            v0=v0,
            r0=r0,
            x0=x0,
            y0=y0,
            Ψ0=Ψ,
            ρ=ρ,
            algorithm=algorithm,
            reltol=reltol,
            abstol=abstol,
        )

        # get finish index
        target_Ψ_rad = Ψ0 + target_Ψ_rad_deviation
        if target_δ_rad < 0
            target_Ψ_rad = Ψ0 - target_Ψ_rad_deviation
        end
        over_index = findfirst(e -> e > target_Ψ_rad, Ψ_list)
        if target_δ_rad < 0
            over_index = findfirst(e -> e < target_Ψ_rad, Ψ_list)
        end
        next_stage_index = length(time_list)
        if isnothing(over_index)
            final_δ_list[start_index:next_stage_index] = δ_list
            final_u_list[start_index:next_stage_index] = u
            final_v_list[start_index:next_stage_index] = v
            final_r_list[start_index:next_stage_index] = r
            final_x_list[start_index:next_stage_index] = x
            final_y_list[start_index:next_stage_index] = y
            final_Ψ_list[start_index:next_stage_index] = Ψ_list
        else
            Ψ = Ψ_list[over_index]
            next_stage_index = over_index + start_index - 1
            final_δ_list[start_index:next_stage_index-1] = δ_list[begin:over_index-1]
            final_u_list[start_index:next_stage_index-1] = u[begin:over_index-1]
            final_v_list[start_index:next_stage_index-1] = v[begin:over_index-1]
            final_r_list[start_index:next_stage_index-1] = r[begin:over_index-1]
            final_x_list[start_index:next_stage_index-1] = x[begin:over_index-1]
            final_y_list[start_index:next_stage_index-1] = y[begin:over_index-1]
            final_Ψ_list[start_index:next_stage_index-1] = Ψ_list[begin:over_index-1]
        end
    end
    final_u_list,
    final_v_list,
    final_r_list,
    final_x_list,
    final_y_list,
    final_Ψ_list,
    final_δ_list
end

function create_model_for_mcmc_sample_mmg(
    data::ShipData,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2;
    ρ=1025.0,
    σ_u_prior_dist=Uniform(0.00, 0.20),
    σ_v_prior_dist=Uniform(0.00, 0.20),
    σ_r_prior_dist=Uniform(0.00, 0.20),
    R_0_dash_prior_dist=Uniform(0.000, 0.100),
    X_vv_dash_prior_dist=Uniform(-0.200, 0.200),
    X_vr_dash_prior_dist=Uniform(-0.223, 0.177),
    X_rr_dash_prior_dist=Uniform(-0.088, 0.032),
    X_vvvv_dash_prior_dist=Uniform(-1.400, 1.400),
    Y_v_dash_prior_dist=Uniform(-0.500, 0.000),
    Y_r_dash_prior_dist=Uniform(-0.100, 0.200),
    Y_vvv_dash_prior_dist=Uniform(-6.000, 2.000),
    Y_vvr_dash_prior_dist=Uniform(-2.500, 1.000),
    Y_vrr_dash_prior_dist=Uniform(-1.500, 0.000),
    Y_rrr_dash_prior_dist=Uniform(-0.120, 0.040),
    N_v_dash_prior_dist=Uniform(-0.200, 0.000),
    N_r_dash_prior_dist=Uniform(-0.100, 0.000),
    N_vvv_dash_prior_dist=Uniform(-0.500, 0.400),
    N_vvr_dash_prior_dist=Uniform(-1.000, 0.000),
    N_vrr_dash_prior_dist=Uniform(-0.300, 0.300),
    N_rrr_dash_prior_dist=Uniform(-0.060, 0.000),
    solver=Tsit5(),
    abstol=1e-6,
    reltol=1e-3
)
    time_obs = data.time
    u_obs = data.u
    v_obs = data.v
    r_obs = data.r
    δ_obs = data.δ
    n_p_obs = data.n_p

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
    x_R,
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
    spl_n_p = Spline1D(time_obs, n_p_obs)
    function MMG!(dX, X, p, t)
        u, v, r, δ, n_p = X
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
        if n_p == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (n_p * D_p)
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
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * n_p^2 * D_p^4 / pi)^2)
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
        X_P = (1.0 - t_P) * ρ * K_T * n_p^2 * D_p^4
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
        N_R = -(x_R + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                    ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dδ = derivative(spl_δ, t)
        dX[5] = dn_p = derivative(spl_n_p, t)
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

    p = (
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
    )

    u0 = 2.29 * 0.512
    v0 = 0.0
    r0 = 0.0
    X0 = [u_obs[1]; v_obs[1]; r_obs[1]; δ_obs[1]; n_p_obs[1]]
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

        p = (
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
        prob = remake(prob1, p=p)
        sol = solve(prob, solver, abstol=abstol, reltol=reltol)
        predicted = sol(time_obs)
        for i in eachindex(predicted)
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
    multivariate_prior_dist;
    ρ=1025.0,
    σ_u_prior_dist=Uniform(0.00, 0.20),
    σ_v_prior_dist=Uniform(0.00, 0.20),
    σ_r_prior_dist=Uniform(0.00, 0.20),
    solver=Tsit5(),
    abstol=1e-6,
    reltol=1e-3
)
    time_obs = data.time
    u_obs = data.u
    v_obs = data.v
    r_obs = data.r
    δ_obs = data.δ
    n_p_obs = data.n_p

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
    x_R,
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
    spl_n_p = Spline1D(time_obs, n_p_obs)
    function MMG!(dX, X, p, t)
        u, v, r, δ, n_p = X
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
        if n_p == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (n_p * D_p)
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
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * n_p^2 * D_p^4 / pi)^2)
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
        X_P = (1.0 - t_P) * ρ * K_T * n_p^2 * D_p^4
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
        N_R = -(x_R + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                    ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dδ = derivative(spl_δ, t)
        dX[5] = dn_p = derivative(spl_n_p, t)
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

    p = (
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
    )

    u0 = 2.29 * 0.512
    v0 = 0.0
    r0 = 0.0
    X0 = [u_obs[1]; v_obs[1]; r_obs[1]; δ_obs[1]; n_p_obs[1]]
    prob1 = ODEProblem(MMG!, X0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    @model function fitMMG(time_obs, obs, prob1)
        σ_u ~ σ_u_prior_dist
        σ_v ~ σ_v_prior_dist
        σ_r ~ σ_r_prior_dist
        p ~ multivariate_prior_dist;
        prob = remake(prob1, p=p)
        sol = solve(prob, solver, abstol=abstol, reltol=reltol)
        predicted = sol(time_obs)
        for i in eachindex(predicted)
            obs[1][i] ~ Normal(predicted[i][1], σ_u) # u
            obs[2][i] ~ Normal(predicted[i][2], σ_v) # v
            obs[3][i] ~ Normal(predicted[i][3], σ_r) # r
        end
    end

    return fitMMG(time_obs, [u_obs, v_obs, r_obs], prob1)
end

function mpc_for_maneuvering_variables(
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams,
    ref_data::BowSternCoords;
    u0=0.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    ρ=1025.0,
    T_all=100.0,
    T_step=1.0,
    Np=1,
    L_bow=3.0,
    L_stern=3.0,
    δ0=1e-5,
    δ_step=5 / 180 * π,
    δ_max=40 / 180 * π,
    δ_min=-40 / 180 * π,
    n_p0=1e-5,
    n_p_step=5.0,
    n_p_max=30.0,
    n_p_min=0.0,
    Q=[1e-3, 1e-3]
    )

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
    x_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    @unpack k_0,
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
    N_rrr_dash = maneuvering_params

    function MMG!(dX, X, δ, n_p)
        u, v, r, x, y, Ψ = X

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
        if n_p == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (n_p * D_p)
        end
        K_T = k_0 + k_1 * J + k_2 * J^2
        β_R = β - l_R * r_dash
        γ_R = γ_R_minus

        γ_R = 0.5 * (γ_R_plus + γ_R_minus) + 0.5 * (γ_R_plus - γ_R_minus) * tanh(100.0 * β_R)

        v_R = U * γ_R * β_R

        u_R = 0.0
        if J == 0.0
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * n_p^2 * D_p^4 / pi)^2)
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
        X_P = (1.0 - t_P) * ρ * K_T * n_p^2 * D_p^4
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
        N_R = -(x_R + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R) * x_G * m +
                    ((Y_H + Y_R) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dx = u * cos(Ψ) - v * sin(Ψ)
        dX[5] = dy = u * sin(Ψ) + v * cos(Ψ)
        dX[6] = dΨ = r
    end

    function rk4_step(X_in, δ, n_p)
        k1 = Vector{Any}(undef, 6)
        MMG!(k1, X_in, δ, n_p)
        X_k2 = X_in + 0.5 * T_step * k1
        k2 = Vector{Any}(undef, 6)
        MMG!(k2, X_k2, δ, n_p)
        X_k3 = X_in + 0.5 * T_step * k2
        k3 = Vector{Any}(undef, 6)
        MMG!(k3, X_k3, δ, n_p)
        X_k4 = X_in + T_step * k3
        k4 = Vector{Any}(undef, 6)
        MMG!(k4, X_k4, δ, n_p)
        X_out = X_in + (T_step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        return X_out
    end
    function discrete_u(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[1]
    end
    function discrete_v(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[2]
    end
    function discrete_r(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[3]
    end
    function discrete_x(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[4]
    end
    function discrete_y(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[5]
    end
    function discrete_Ψ(u, v, r, x, y, Ψ, δ, n_p)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p)
        return X_out[6]
    end
    
    total_steps = Int(T_all / T_step)
    X0_mpc = [u0, v0, r0, x0, y0, Ψ0, x0 + L_bow * cos(Ψ0), y0 + L_bow * sin(Ψ0), x0 + L_stern * cos(Ψ0 + pi), y0 + L_stern * sin(Ψ0 + pi)]
    U0_mpc = [δ0, n_p0]
    
    results_list = []

    @showprogress for k in 1:total_steps
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        register(model, :discrete_u, 8, discrete_u; autodiff = true)
        register(model, :discrete_v, 8, discrete_v; autodiff = true)
        register(model, :discrete_r, 8, discrete_r; autodiff = true)
        register(model, :discrete_x, 8, discrete_x; autodiff = true)
        register(model, :discrete_y, 8, discrete_y; autodiff = true)
        register(model, :discrete_Ψ, 8, discrete_Ψ; autodiff = true)

        @variable(model, X_mpc[1:Np+1, 1:10])
        @variable(model, U_mpc[1:Np, 1:2])
        @constraint(model, X_mpc[1, 1:10] .== X0_mpc[1:10])
        @constraint(model, δ_min .<= U_mpc[1:Np, 1] .<= δ_max)
        @constraint(model, n_p_min .<= U_mpc[1:Np, 2] .<= n_p_max)        
        @constraint(model, -δ_step .<= U_mpc[1, 1] - U0_mpc[1] .<= δ_step)
        @constraint(model, -n_p_step .<= U_mpc[1, 2] - U0_mpc[2] .<= n_p_step)
        for i in 1:Np-1
            @constraint(model, -δ_step .<= U_mpc[i+1, 1] - U_mpc[i, 1] .<= δ_step)
            @constraint(model, -n_p_step .<= U_mpc[i+1, 2] - U_mpc[i, 2] .<= n_p_step)
        end

        for i in 1:Np+1
            set_start_value.(X_mpc[i, 1:10], X0_mpc[1:10])
        end
        for i in 1:Np
            set_start_value.(U_mpc[i, 1:2], U0_mpc[1:2])
        end  

        for i in 1:Np
            @NLconstraint(model, X_mpc[i+1, 1] == discrete_u(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 2] == discrete_v(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 3] == discrete_r(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 4] == discrete_x(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 5] == discrete_y(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 6] == discrete_Ψ(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], U_mpc[i, 1], U_mpc[i, 2]))
            @NLconstraint(model, X_mpc[i+1, 7] == X_mpc[i+1, 4] + L_bow * cos(X_mpc[i+1, 6]))
            @NLconstraint(model, X_mpc[i+1, 8] == X_mpc[i+1, 5] + L_bow * sin(X_mpc[i+1, 6]))
            @NLconstraint(model, X_mpc[i+1, 9] == X_mpc[i+1, 4] + L_stern * cos(X_mpc[i+1, 6] + pi))
            @NLconstraint(model, X_mpc[i+1, 10] == X_mpc[i+1, 5] + L_stern * sin(X_mpc[i+1, 6]+ pi))
        end

        @NLobjective(model, Min, 
            sum((X_mpc[i, 7]-ref_data.x1[k+i-1])^2 + (X_mpc[i, 8]-ref_data.y1[k+i-1])^2 + (X_mpc[i, 9]-ref_data.x2[k+i-1])^2 + (X_mpc[i, 10]-ref_data.y2[k+i-1])^2  for i in 2:Np+1) 
            + sum(Q[1] * (U_mpc[i, 1]-U_mpc[i-1, 1])^2 + Q[2] * (U_mpc[i, 2]-U_mpc[i-1, 2])^2 for i in 2:Np))
            + (Q[1] * (U_mpc[1, 1]-U0_mpc[1])^2 + Q[2] * (U_mpc[1, 2]-U0_mpc[2])^2)
        
        optimize!(model)
        
        U0_mpc = JuMP.value.(U_mpc[1, 1:2])
        push!(results_list, (
            time = (k - 1) * T_step,
            state = X0_mpc,
            control = U0_mpc
        ))
        X0_mpc = JuMP.value.(X_mpc[2, 1:10])
    end

    push!(results_list, (
        time = total_steps * T_step,
        state = X0_mpc,
        control = [NaN, NaN]
    ))

    time = [res.time for res in results_list]
    u = [res.state[1] for res in results_list]
    v = [res.state[2] for res in results_list]
    r = [res.state[3] for res in results_list]
    x = [res.state[4] for res in results_list]
    y = [res.state[5] for res in results_list]
    Ψ = [res.state[6] for res in results_list]
    x1 = [res.state[7] for res in results_list]
    y1 = [res.state[8] for res in results_list]
    x2 = [res.state[9] for res in results_list]
    y2 = [res.state[10] for res in results_list]
    δ = [res.control[1] for res in results_list]
    n_p = [res.control[2] for res in results_list]

    return time, u, v, r, x, y, Ψ, x1, y1, x2, y2, δ, n_p
end

function mpc_for_external_force(
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams,
    ref_data::BowSternCoords,
    δ_ref::Vector{Float64},
    n_p_ref::Vector{Float64};
    u0=0.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    ρ=1025.0,
    T_all=100.0,
    T_step=1.0,
    Np=1,
    L_bow=3.0,
    L_stern=3.0,
    X_F0=1e-5,
    X_F_max=1e2,
    X_F_min=-1e2,
    Y_F0=1e-5,
    Y_F_max=1e2,
    Y_F_min=-1e2,
    N_F0=1e-5,
    N_F_max=1e2,
    N_F_min=-1e2,
    Q=[1e-6, 1e-6, 1e-6]
    )

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
    x_R,
    a_H,
    x_H,
    γ_R_minus,
    γ_R_plus,
    l_R,
    κ,
    t_P,
    w_P0,
    x_P = basic_params

    @unpack k_0,
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
    N_rrr_dash = maneuvering_params

    function MMG!(dX, X, δ, n_p, X_F, Y_F, N_F)
        u, v, r, x, y, Ψ = X

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
        if n_p == 0.0
            J = 0.0
        else
            J = (1 - w_P) * u / (n_p * D_p)
        end
        K_T = k_0 + k_1 * J + k_2 * J^2
        β_R = β - l_R * r_dash
        γ_R = γ_R_minus

        γ_R = 0.5 * (γ_R_plus + γ_R_minus) + 0.5 * (γ_R_plus - γ_R_minus) * tanh(100.0 * β_R)

        v_R = U * γ_R * β_R

        u_R = 0.0
        if J == 0.0
            u_R = sqrt(η * (κ * ϵ * 8.0 * k_0 * n_p^2 * D_p^4 / pi)^2)
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
        X_P = (1.0 - t_P) * ρ * K_T * n_p^2 * D_p^4
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
        N_R = -(x_R + a_H * x_H) * F_N * cos(δ)
        dX[1] = du = ((X_H + X_R + X_P + X_F) + (m + m_y) * v * r + x_G * m * (r^2)) / (m + m_x)
        dX[2] =
            dv =
                (
                    (x_G^2) * (m^2) * u * r - (N_H + N_R + N_F) * x_G * m +
                    ((Y_H + Y_R + Y_F) - (m + m_x) * u * r) * (I_zG + J_z + (x_G^2) * m)
                ) / ((I_zG + J_z + (x_G^2) * m) * (m + m_y) - (x_G^2) * (m^2))
        dX[3] = dr = (N_H + N_R + N_F - x_G * m * (dv + u * r)) / (I_zG + J_z + (x_G^2) * m)
        dX[4] = dx = u * cos(Ψ) - v * sin(Ψ)
        dX[5] = dy = u * sin(Ψ) + v * cos(Ψ)
        dX[6] = dΨ = r
    end

    function rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        k1 = Vector{Any}(undef, 6)
        MMG!(k1, X_in, δ, n_p, X_F, Y_F, N_F)
        X_k2 = X_in + 0.5 * T_step * k1
        k2 = Vector{Any}(undef, 6)
        MMG!(k2, X_k2, δ, n_p, X_F, Y_F, N_F)
        X_k3 = X_in + 0.5 * T_step * k2
        k3 = Vector{Any}(undef, 6)
        MMG!(k3, X_k3, δ, n_p, X_F, Y_F, N_F)
        X_k4 = X_in + T_step * k3
        k4 = Vector{Any}(undef, 6)
        MMG!(k4, X_k4, δ, n_p, X_F, Y_F, N_F)
        X_out = X_in + (T_step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        return X_out
    end
    function discrete_u(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[1]
    end
    function discrete_v(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[2]
    end
    function discrete_r(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[3]
    end
    function discrete_x(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[4]
    end
    function discrete_y(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[5]
    end
    function discrete_Ψ(u, v, r, x, y, Ψ, δ, n_p, X_F, Y_F, N_F)
        X_in = [u, v, r, x, y, Ψ]
        X_out = rk4_step(X_in, δ, n_p, X_F, Y_F, N_F)
        return X_out[6]
    end

    total_steps = Int(T_all / T_step)
    X0_mpc = [u0, v0, r0, x0, y0, Ψ0, x0 + L_bow * cos(Ψ0), y0 + L_bow * sin(Ψ0), x0 + L_stern * cos(Ψ0 + pi), y0 + L_stern * sin(Ψ0 + pi)]
    U0_mpc = [X_F0, Y_F0, N_F0]

    results_list = []

    @showprogress for k in 1:total_steps
        model = Model(Ipopt.Optimizer)
        set_silent(model)
        register(model, :discrete_u, 11, discrete_u; autodiff = true)
        register(model, :discrete_v, 11, discrete_v; autodiff = true)
        register(model, :discrete_r, 11, discrete_r; autodiff = true)
        register(model, :discrete_x, 11, discrete_x; autodiff = true)
        register(model, :discrete_y, 11, discrete_y; autodiff = true)
        register(model, :discrete_Ψ, 11, discrete_Ψ; autodiff = true)

        @variable(model, X_mpc[1:Np+1, 1:10])
        @variable(model, U_mpc[1:Np, 1:3])
        @constraint(model, X_mpc[1, 1:10] .== X0_mpc[1:10])
        @constraint(model, X_F_min .<= U_mpc[1:Np, 1] .<= X_F_max)
        @constraint(model, Y_F_min .<= U_mpc[1:Np, 2] .<= Y_F_max)
        @constraint(model, N_F_min .<= U_mpc[1:Np, 3] .<= N_F_max)

        for i in 1:Np+1
            set_start_value.(X_mpc[i, 1:10], X0_mpc[1:10])
        end
        for i in 1:Np
            set_start_value.(U_mpc[i, 1:3], U0_mpc[1:3])
        end  

        for i in 1:Np
            @NLconstraint(model, X_mpc[i+1, 1] == discrete_u(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 2] == discrete_v(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 3] == discrete_r(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 4] == discrete_x(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 5] == discrete_y(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 6] == discrete_Ψ(X_mpc[i, 1], X_mpc[i, 2], X_mpc[i, 3], X_mpc[i, 4], X_mpc[i, 5], X_mpc[i, 6], δ_ref[k+i-1], n_p_ref[k+i-1], U_mpc[i, 1], U_mpc[i, 2], U_mpc[i, 3]))
            @NLconstraint(model, X_mpc[i+1, 7] == X_mpc[i+1, 4] + L_bow * cos(X_mpc[i+1, 6]))
            @NLconstraint(model, X_mpc[i+1, 8] == X_mpc[i+1, 5] + L_bow * sin(X_mpc[i+1, 6]))
            @NLconstraint(model, X_mpc[i+1, 9] == X_mpc[i+1, 4] + L_stern * cos(X_mpc[i+1, 6] + pi))
            @NLconstraint(model, X_mpc[i+1, 10] == X_mpc[i+1, 5] + L_stern * sin(X_mpc[i+1, 6]+ pi))
        end

        @NLobjective(model, Min, 
            sum((X_mpc[i, 7]-ref_data.x1[k+i-1])^2 + (X_mpc[i, 8]-ref_data.y1[k+i-1])^2 + (X_mpc[i, 9]-ref_data.x2[k+i-1])^2 + (X_mpc[i, 10]-ref_data.y2[k+i-1])^2  for i in 2:Np+1)  
            + sum(Q[1]*(U_mpc[i, 1]-U_mpc[i-1, 1])^2 + Q[2]*(U_mpc[i, 2]-U_mpc[i-1, 2])^2 + Q[3]*(U_mpc[i, 3]-U_mpc[i-1, 3])^2 for i in 2:Np))
            + (Q[1]*(U_mpc[1, 1]-U0_mpc[1])^2 + Q[2]*(U_mpc[1, 2]-U0_mpc[2])^2 + Q[3]*(U_mpc[1, 3]-U0_mpc[3])^2)
        
        optimize!(model)

        U0_mpc = JuMP.value.(U_mpc[1, 1:3])
        push!(results_list, (
            time = (k - 1) * T_step,
            state = [X0_mpc..., δ_ref[k], n_p_ref[k]],
            control = U0_mpc
        ))
        X0_mpc = JuMP.value.(X_mpc[2, 1:10])
    end

    push!(results_list, (
        time = total_steps * T_step,
        state = [X0_mpc..., δ_ref[total_steps+1], n_p_ref[total_steps+1]],
        control = [NaN, NaN, NaN]
    ))

    time = [res.time for res in results_list]
    u = [res.state[1] for res in results_list]
    v = [res.state[2] for res in results_list]
    r = [res.state[3] for res in results_list]
    x = [res.state[4] for res in results_list]
    y = [res.state[5] for res in results_list]
    Ψ = [res.state[6] for res in results_list]
    x1 = [res.state[7] for res in results_list]
    y1 = [res.state[8] for res in results_list]
    x2 = [res.state[9] for res in results_list]
    y2 = [res.state[10] for res in results_list]
    δ = [res.state[11] for res in results_list]
    n_p = [res.state[12] for res in results_list]
    X_F = [res.control[1] for res in results_list]
    Y_F = [res.control[2] for res in results_list]
    N_F = [res.control[3] for res in results_list]
    
    return time, u, v, r, x, y, Ψ, x1, y1, x2, y2, δ, n_p, X_F, Y_F, N_F
end