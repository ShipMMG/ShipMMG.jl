module ShipMMG

    using DifferentialEquations
    using ParameterizedFunctions
    using Dierckx
    using LinearAlgebra
    using Sundials

    include("kt.jl")
    export kt_simulate
    
end
