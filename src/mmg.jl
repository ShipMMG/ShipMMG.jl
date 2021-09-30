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

function mmg_3dof_simulate(
    time_list,
    npm_list,
    δ_list,
    basic_params::Mmg3DofBasicParams,
    maneuvering_params::Mmg3DofManeuveringParams;
    u0::Float64 = 0.0,
    v0::Float64 = 0.0,
    r0::Float64 = 0.0,
    ρ::Float64 = 1.025,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)
    simulate(
        time_list,
        npm_list,
        δ_list,
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
        u0 = u0,
        v0 = v0,
        r0 = r0,
        ρ = ρ,
        algorithm = algorithm,
        reltol = reltol,
        abstol = abstol,
    )
end

function simulate(
    time_list,
    npm_list,
    δ_list,
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
    N_rrr_dash::Float64;
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

    function mmg_3dof_eom!(dX, X, maneuvering_params_list, t)
        u, v, r, δ, npm = X
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
        N_rrr_dash = maneuvering_params_list

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

    X0 = [u0; v0; r0; δ_list[1]; npm_list[1]]
    maneuvering_params_list = [
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
    ]
    prob = ODEProblem(
        mmg_3dof_eom!,
        X0,
        (time_list[1], time_list[end]),
        maneuvering_params_list,
    )
    sol = solve(prob, algorithm, reltol = reltol, abstol = abstol)
    results = hcat(sol.u...)
    time = sol.t
    u = results[1, :]
    v = results[2, :]
    r = results[3, :]
    δ = results[4, :]
    npm = results[5, :]
    time, u, v, r, δ, npm
end
