module ShipMMG

    using DifferentialEquations
    using ParameterizedFunctions
    using Dierckx
    using LinearAlgebra
    using Sundials
    using Plots
    using LaTeXStrings

    include("kt.jl")
    export kt_simulate, draw_gif_result, rotate_pos, square

    function draw_gif_result(time, x, y, ψ, shape, file_path, fps=10)
        anim = @animate for i in 1:length(time)
            plot(x, y, label="", xlabel=L"x \textrm{[m]}", ylabel=L"y \textrm{[m]}", linestyle = :dot, aspect_ratio=:equal)
            ship = square(x[i], y[i], shape, ψ[i])
            plot!(Shape(ship[1], ship[2]), label="")
            scatter!([x[i]], [y[i]], seriestype = :scatter, title="time = $(time[i])", label="")
        end
        gif(anim, file_path, fps=fps)
    end

    function rotate_pos(pos, angle)
        rotate_matrix = [cos(angle) -sin(angle) ; sin(angle) cos(angle)]
        rotate_matrix * pos
    end

    function square(center_x, center_y, shape, angle)
        square_xy = [
            rotate_pos([shape[1], shape[2]]/2,angle)+[center_x,center_y],
            rotate_pos([-shape[1], shape[2]]/2,angle)+[center_x,center_y],
            rotate_pos([-shape[1], -shape[2]]/2,angle)+[center_x,center_y],
            rotate_pos([shape[1], -shape[2]]/2,angle)+[center_x,center_y],
        ]
        xy = hcat(square_xy...)
        [xy[1,:], xy[2,:]]
    end
    
end
