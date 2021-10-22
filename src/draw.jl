function calc_position(
    time_vec::Vector{Float64},
    u_vec::Vector{Float64},
    v_vec::Vector{Float64},
    r_vec::Vector{Float64};
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

"""
    draw_gif_result(time, x, y, ψ, shape, file_path, [, fps]) -> gif

Draw the gif animation from simulation result.
"""
function draw_gif_result(time, x, y, ψ, shape, file_path, fps = 10)
    anim = @animate for i = 1:length(time)
        plot(
            x,
            y,
            label = "",
            xlabel = L"x \textrm{[m]}",
            ylabel = L"y \textrm{[m]}",
            linestyle = :dot,
            aspect_ratio = :equal,
        )
        ship = square(x[i], y[i], shape, ψ[i])
        plot!(Shape(ship[1], ship[2]), label = "")
        scatter!(
            [x[i]],
            [y[i]],
            seriestype = :scatter,
            title = "time = $(time[i])",
            label = "",
        )
    end
    gif(anim, file_path, fps = fps)
end

function rotate_pos(pos, angle)
    rotate_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    rotate_matrix * pos
end

function square(center_x, center_y, shape, angle)
    square_xy = [
        rotate_pos([shape[1], shape[2]] / 2, angle) + [center_x, center_y],
        rotate_pos([-shape[1], shape[2]] / 2, angle) + [center_x, center_y],
        rotate_pos([-shape[1], -shape[2]] / 2, angle) + [center_x, center_y],
        rotate_pos([shape[1], -shape[2]] / 2, angle) + [center_x, center_y],
    ]
    xy = hcat(square_xy...)
    [xy[1, :], xy[2, :]]
end