module ShipMMG

using DifferentialEquations
using ParameterizedFunctions
using Dierckx
using LinearAlgebra
using Sundials
using Plots
using LaTeXStrings

include("kt.jl")
export kt_simulate

include("mmg.jl")
export Mmg3DofBasicParams, Mmg3DofManeuveringParams, mmg_3dof_simulate

include("draw.jl")
export draw_gif_result, calc_position

end
