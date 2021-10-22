"""
    kt_model!(dX, X, p, t)

KT model on DifferentialEquations.ODEProblem. Update `dX`.

# Arguments
- `dX`: [du, dr, dδ]
- `X`: the initial state values. [`u`, `r`, `δ`].
- `p`: the parameters and δ spline info. [`K`, `T`, `spl_δ`].
- `t`: the time.

# Examples

```julia-rep1
julia> K_log = 0.155  # [1/s]
julia> T_log = 80.5  # [s]
julia> u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
julia> duration = 50  # [s]
julia> sampling = 1001
julia> time_list = range(0.0,stop=duration,length=sampling)
julia> Ts = 50.0
julia> δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
julia> spl_δ = Spline1D(time_list, δ_list)
julia> X0 = [u0; 0.0; δ_list[1]]
julia> p = [K_log, T_log, spl_δ]
julia> prob = ODEProblem(kt_model!, X0, (time_list[1], time_list[end]), p)
julia> sol = solve(prob, Tsit5(), saveat=time_list[2] - time_list[1])
```
"""
function kt_model!(dX, X, p, t)
    u, r, δ = X
    K, T, spl_δ = p
    dX[1] = du = 0.0 # du = 
    dX[2] = dr = 1.0 / T * (-r + K * δ) # dr = 
    dX[3] = dδ = derivative(spl_δ, t) # dδ = 
end

"""
    kt_simulate(time_list, δ_list, K, T, u0, [, r0, algorithm, reltol, abstol])

Returns the KT simulation results.

# Arguments
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `K::Float64`: the K Parameter.
- `T::Float64`: the T Parameter.
- `u0::Float64`: the constant surge velocity.
- `r0::Float64=0.0`: the inital rate of turn [rad/s].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KT simulation.

```julia-rep1
julia> K_log = 0.155  # [1/s]
julia> T_log = 80.5  # [s]
julia> u0 = 20 * (1852.0 / 3600)  # [m/s] (knot * 1852/3600)
julia> duration = 50  # [s]
julia> sampling = 1001
julia> time_list = range(0.0,stop=duration,length=sampling)
julia> Ts = 50.0
julia> δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
julia> results=kt_simulate(time_list, δ_list, K_log, T_log, u0)
```
"""
function kt_simulate(
    time_list,
    δ_list,
    K::Float64,
    T::Float64,
    u0::Float64,
    r0::Float64 = 0.0,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)

    spl_δ = Spline1D(time_list, δ_list)

    X0 = [u0; r0; δ_list[1]]
    p = [K, T, spl_δ]
    prob = ODEProblem(kt_model!, X0, (time_list[1], time_list[end]), p)
    sol = solve(
        prob,
        algorithm,
        reltol = reltol,
        abstol = abstol,
        saveat = time_list[2] - time_list[1],
    )

    results = hcat(sol.u...)
    time = sol.t
    u = results[1, :]
    r = results[2, :]
    δ = results[3, :]
    time, u, r, δ
end
