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
    r_dot = f(spl_r, time_vec)
    A = hcat(δ_vec, r_vec)
    Θ = A \ r_dot
    T = -1.0 / Θ[2]
    K = Θ[1] * T
    K, T
end

function estimate_mmg_approx_lsm(
    data::ShipData,
    basic_params::Mmg3DofBasicParams,
    k_0,
    k_1,
    k_2;
    ρ = 1.025,
)
    @unpack time_vec, u, v, r, x, y, ψ, δ, npm = data
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
    for i = 1:length(time_vec)
        if isapprox(U[i], 0.0) == false
            β[i] = asin(-(v[i] - r[i] * x_G) / U[i])
        end
    end

    v_dash = U .^ 0
    for i = 1:length(time_vec)
        if isapprox(U[i], 0.0) == false
            v_dash[i] = v[i] / U[i]
        end
    end

    r_dash = U .^ 0
    for i = 1:length(time_vec)
        if isapprox(U[i], 0.0) == false
            r_dash[i] = r[i] * L_pp / U[i]
        end
    end

    w_P = w_P0 * exp.(.-4.0 * (β - r_dash .* x_P) .^ 2)

    J = U .^ 0
    for i = 1:length(time_vec)
        if isapprox(npm[i], 0.0) == false
            J[i] = (1.0 - w_P[i]) * u[i] / (npm[i] * D_p)
        end
    end

    K_T = J .^ 2 .* k_2 + J .* k_1 .+ k_0
    β_R = β - r_dash .* l_R

    γ_R = U .^ 0
    for i = 1:length(time_vec)
        if β_R[i] < 0.0
            γ_R[i] = γ_R_minus
        else
            γ_R[i] = γ_R_plus
        end
    end

    v_R = U .* γ_R .* β_R

    u_R = U .^ 0
    for i = 1:length(time_vec)
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
    spl_u = Spline1D(time_vec, u)
    spl_v = Spline1D(time_vec, v)
    spl_r = Spline1D(time_vec, r)
    f = (spl, x) -> derivative(spl, x)
    du = f(spl_u, time_vec)
    dv = f(spl_v, time_vec)
    dr = f(spl_r, time_vec)


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
