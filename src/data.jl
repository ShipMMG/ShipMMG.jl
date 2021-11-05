@with_kw struct ShipData{T}
    time::T
    u::T
    v::T
    r::T
    x::T
    y::T
    ψ::T
    δ::T
    npm::T
end

function calc_position(
    time_vec,
    u_vec,
    v_vec,
    r_vec;
    x0::Float64 = 0.0,
    y0::Float64 = 0.0,
    ψ0::Float64 = 0.0,
)
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