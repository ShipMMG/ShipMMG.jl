module ShipMMG

using DifferentialEquations
using ParameterizedFunctions
using Dierckx
using Parameters
using Distributions
using Turing
using ForwardDiff
using LinearAlgebra
using JuMP
using Ipopt

include("data.jl")
export calc_position,
    ShipData,
    EnvironmentData,
    PointPairCoords,
    get_KVLCC2_L7_basic_params,
    get_KVLCC2_L7_maneuvering_params,
    get_KVLCC2_L7_params,
    get_example_ship_wind_force_moment_params,
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
    mmg_3dof_simulate,
    mmg_3dof_model!,
    mmg_3dof_zigzag_test,
    estimate_mmg_approx_lsm,
    estimate_mmg_approx_lsm_time_window_sampling,
    create_model_for_mcmc_sample_mmg,
    mpc_for_maneuvering_variables,
    mpc_for_external_force

include("mmg_wind.jl")
export Mmg3DofBasicParams,
    Mmg3DofManeuveringParams,
    Mmg3DofWindForceMomentParams,
    apparent_wind_speed_and_angle,
    mmg_3dof_simulate,
    mmg_3dof_wind_model!,
    mmg_3dof_zigzag_test,
    estimate_mmg_approx_lsm,
    estimate_mmg_approx_lsm_time_window_sampling,
    create_model_for_mcmc_sample_mmg,
    wind_force_and_moment_coefficients

end
