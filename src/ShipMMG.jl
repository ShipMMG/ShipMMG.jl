module ShipMMG

using DifferentialEquations
using ParameterizedFunctions
using Dierckx
using Plots
using Parameters

include("data.jl")
export calc_position, ShipData

include("kt.jl")
export kt_simulate, kt_model!, kt_zigzag_test

include("mmg.jl")
export Mmg3DofBasicParams,
    Mmg3DofManeuveringParams, mmg_3dof_simulate, mmg_3dof_model!, mmg_3dof_zigzag_test

include("draw.jl")
export draw_gif_result, calc_position

end
