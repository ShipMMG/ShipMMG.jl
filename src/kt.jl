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
    kt_simulate(time_list, δ_list, K, T, [, u0, r0, algorithm, reltol, abstol])

Returns the KT simulation results.

# Arguments
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `K::Float64`: the K Parameter.
- `T::Float64`: the T Parameter.
- `u0::Float64=0.0`: the constant surge velocity.
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
julia> time, u, r, δ = kt_simulate(time_list, δ_list, K_log, T_log, u0=u0)
```
"""
function kt_simulate(
    time_list,
    δ_list,
    K::Float64,
    T::Float64;
    u0::Float64 = 0.0,
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

"""
    kt_zigzag_test(K, T, target_δ_rad, target_ψ_rad_deviation, time_second_interval, end_time_second, [, u0, r0, ψ0, δ0, δ_rad_rate, start_time_second, algorithm, reltol, abstol])

Returns the KT simulation results.

# Arguments
- `K::Float64`: the K Parameter.
- `T::Float64`: the T Parameter.
- `target_δ_rad::Float64`: target rudder angle of zigzag test.
- `target_ψ_rad_deviation::Float64`: target azimuth deviation of zigzag test.
- `time_second_interval::Float64`: time interval of output result.
- `end_time_second::Float64`: the end time of simulation.
- `u0::Float64=0.0`: the initial & constant surge velocity.
- `r0::Float64=0.0`: the initial rate of turn [rad/s].
- `ψ0::Float64=0.0`: the initial azimuth.
- `δ0::Float64=0.0`: the initial rudder angle.
- `δ_rad_rate::Float64=10.0*π/180`: the change rate of rudder angle [rad/s]. 
- `start_time_second::Float64=0.0`: the start time of simulation.
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KT simulation.

```julia-rep1
julia> K = 0.155
julia> T = 80.5
julia> target_δ_rad = 30.0 * π / 180.0
julia> target_ψ_rad_deviation = 10.0 * π / 180.0
julia> time_second_interval = 0.01
julia> end_time_second = 500.0
julia> u0 = 20 * (1852.0 / 3600) # [knot]
julia> δ_list, r_list, ψ_list = kt_zigzag_test(
        K,
        T,
        target_δ_rad,
        target_ψ_rad_deviation,
        time_second_interval,
        end_time_second,
        u0 = u0,
    )
```
"""
function kt_zigzag_test(
    K::Float64,
    T::Float64,
    target_δ_rad::Float64,
    target_ψ_rad_deviation::Float64,
    time_second_interval::Float64,
    end_time_second::Float64;
    u0::Float64 = 0.0,
    r0::Float64 = 0.0,
    ψ0::Float64 = 0.0,
    δ0::Float64 = 0.0,
    δ_rad_rate::Float64 = 10.0 * π / 180,
    start_time_second::Float64 = 0.0,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)
    target_ψ_rad_deviation = abs(target_ψ_rad_deviation)

    time_list = start_time_second:time_second_interval:end_time_second
    final_δ_list = zeros(length(time_list))
    final_r_list = zeros(length(time_list))
    final_ψ_list = zeros(length(time_list))

    next_stage_index = 1
    target_δ_rad = -target_δ_rad  # for changing in while loop
    ψ = ψ0

    while next_stage_index < length(time_list)
        target_δ_rad = -target_δ_rad
        start_index = next_stage_index

        # Make delta list
        δ_list = zeros(length(time_list) - start_index + 1)
        if start_index == 1
            δ_list[1] = δ0
            r0 = r0
        else
            δ_list[1] = final_δ_list[start_index-1]
            r0 = final_r_list[start_index-1]
        end

        for i = (start_index+1):length(time_list)
            Δt = time_list[i] - time_list[i-1]
            if target_δ_rad > 0
                δ = δ_list[i-start_index] + δ_rad_rate * Δt
                if δ >= target_δ_rad
                    δ = target_δ_rad
                end
                δ_list[i-start_index+1] = δ
            elseif target_δ_rad <= 0
                δ = δ_list[i-start_index] - δ_rad_rate * Δt
                if δ <= target_δ_rad
                    δ = target_δ_rad
                end
                δ_list[i-start_index+1] = δ
            end
        end

        time, u, r, δ =
            kt_simulate(time_list[start_index:end], δ_list, K, T, u0 = u0, r0 = r0)
        v = zeros(length(u))
        x, y, ψ_list = calc_position(time, u, v, r, x0 = 0.0, y0 = 0.0, ψ0 = ψ)
        # get finish index
        target_ψ_rad = ψ0 + target_ψ_rad_deviation
        if target_δ_rad < 0
            target_ψ_rad = ψ0 - target_ψ_rad_deviation
        end
        over_index = findfirst(e -> e >= target_ψ_rad, ψ_list)
        if target_δ_rad < 0
            over_index = findfirst(e -> e <= target_ψ_rad, ψ_list)
        end
        next_stage_index = length(time_list)
        if isnothing(over_index)
            final_δ_list[start_index:next_stage_index] = δ_list
            final_r_list[start_index:next_stage_index] = r
            final_ψ_list[start_index:next_stage_index] = ψ_list
        else
            ψ = ψ_list[over_index]
            next_stage_index = over_index + start_index - 1
            final_δ_list[start_index:next_stage_index] = δ_list[begin:over_index]
            final_r_list[start_index:next_stage_index] = r[begin:over_index]
            final_ψ_list[start_index:next_stage_index] = ψ_list[begin:over_index]
        end
    end
    final_δ_list, final_r_list, final_ψ_list
end

