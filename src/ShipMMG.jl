module ShipMMG

using DifferentialEquations
using ParameterizedFunctions
using Dierckx
using Plots
using Parameters
using Distributions
using Turing
using ForwardDiff

include("data.jl")
export calc_position,
    ShipData,
    get_KVLCC2_L7_basic_params,
    get_KVLCC2_L7_maneuvering_params,
    get_KVLCC2_L7_params,
    nuts_sampling_single_thread,
    nuts_sampling_multi_threads

include("kt.jl")
export kt_simulate,
    kt_model!,
    kt_zigzag_test,
    estimate_kt_lsm,
    estimate_kt_lsm_time_window_sampling,
    create_model_for_mcmc_sample_kt

include("mmg.jl")
export Mmg3DofBasicParams,
    Mmg3DofManeuveringParams,
    Mmg3DofWindForceMomentParams,
    apparent_wind_speed_and_angle,
    mmg_3dof_simulate,
    mmg_3dof_model!,
    mmg_3dof_zigzag_test,
    estimate_mmg_approx_lsm,
    estimate_mmg_approx_lsm_time_window_sampling,
    create_model_for_mcmc_sample_mmg

include("draw.jl")
export draw_gif_result, calc_position

end
