"""
    kt_model!(dX, X, p, t)

KT model on DifferentialEquations.ODEProblem. Update `dX`.

# Arguments
- `dX`: [dr, dx, dy, dΨ, dδ]
- `X`: the initial state values. [`r`, `x`, `y`, `Ψ`, `δ`].
- `p`: the parameters and δ spline info. [`K`, `T`, `spl_δ`].
- `t`: the time.

# Examples

```julia-rep1
julia> K_log = 0.155  # [1/s]
julia> T_log = 80.5  # [s]
julia> duration = 50  # [s]
julia> sampling = 1001
julia> time_list = range(0.0,stop=duration,length=sampling)
julia> Ts = 50.0
julia> δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
julia> spl_δ = Spline1D(time_list, δ_list)
julia> X0 = [10.0; 0.0; 0.0; 0.0; 0.0, 0.0; δ_list[1]]
julia> p = [K_log, T_log, spl_δ]
julia> prob = ODEProblem(kt_model!, X0, (time_list[1], time_list[end]), p)
julia> sol = solve(prob, Tsit5(), saveat=time_list[2] - time_list[1])
```
"""
function kt_model!(dX, X, p, t)
    u, v, r, x, y, Ψ, δ = X
    K, T, spl_δ = p
    dX[1] = du = 0.0
    dX[2] = dv = 0.0
    dX[3] = dr = 1.0 / T * (-r + K * δ)
    dX[4] = dx = u * cos(Ψ)
    dX[5] = dy = u * sin(Ψ)
    dX[6] = dΨ = r
    dX[7] = dδ = derivative(spl_δ, t)
end

"""
    kt_simulate(time_list, δ_list, K, T, [, u0, v0, r0, x0, y0, Ψ0, algorithm, reltol, abstol]) -> u, v, r, x, y, Ψ, δ

Returns the KT simulation results.

# Arguments
- `K`: the K Parameter.
- `T`: the T Parameter.
- `time_list`: the list of simulatino time.
- `δ_list`: the list of rudder angle [rad].
- `u0=10.0`: the initial x (surge) velocity.
- `v0=0.0`: the initial y (sway) velocity.
- `r0=0.0`: the inital rate of turn [rad/s].
- `x0=0.0`: the initial x (surge) position.
- `y0=0.0`: the initial y (sway) position.
- `Ψ0=0.0`: the initial Ψ (yaw) azimuth [rad].
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KT simulation.

```julia-rep1
julia> K_log = 0.155  # [1/s]
julia> T_log = 80.5  # [s]
julia> duration = 50  # [s]
julia> sampling = 1001
julia> time_list = range(0.0,stop=duration,length=sampling)
julia> Ts = 50.0
julia> δ_list = 10.0 * pi / 180.0  * sin.(2.0 * pi / Ts * time_list) # [rad]
julia> u, v, r, x, y, Ψ, δ = kt_simulate(K_log, T_log, time_list, δ_list, u0 = 10.0)
```
"""
function kt_simulate(
    K,
    T,
    time_list,
    δ_list;
    u0=10.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    algorithm=Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
)
    spl_δ = Spline1D(time_list, δ_list)

    X0 = [u0; v0; r0; x0; y0; Ψ0; δ_list[1]]
    p = (K, T, spl_δ)
    prob = ODEProblem(kt_model!, X0, (time_list[1], time_list[end]), p)
    sol = solve(prob, algorithm, reltol=reltol, abstol=abstol)
    sol_timelist = sol(time_list)
    results = hcat(sol_timelist.u...)
    u = results[1, :]
    v = results[2, :]
    r = results[3, :]
    x = results[4, :]
    y = results[5, :]
    Ψ = results[6, :]
    δ = results[7, :]
    u, v, r, x, y, Ψ, δ
end

"""
    kt_zigzag_test(K, T, target_δ_rad, target_Ψ_rad_deviation, time_second_interval, end_time_second, [, u0, v0, r0, x0, y0, Ψ0, δ0, δ_rad_rate, start_time_second, algorithm, reltol, abstol]) -> u, v, r, x, y, Ψ, δ

Returns the KT simulation results.

# Arguments
- `K`: the K Parameter.
- `T`: the T Parameter.
- `time_list`: the list of simulatino time.
- `target_δ_rad`: target rudder angle of zigzag test.
- `target_Ψ_rad_deviation`: target azimuth deviation of zigzag test.
- `u0=10.0`: the initial x (surge) velocity.
- `v0=0.0`: the initial y (sway) velocity.
- `r0=0.0`: the initial rate of turn [rad/s].
- `x0=0.0`: the initial x (surge) position.
- `y0=0.0`: the initial y (sway) position.
- `Ψ0=0.0`: the initial Ψ (yaw) azimuth [rad].
- `δ0=0.0`: the initial rudder angle.
- `δ_rad_rate=10.0*π/180`: the change rate of rudder angle [rad/s]. 
- `algorithm=Tsit5()`: the parameter of DifferentialEquations.ODEProblem.solve()
- `reltol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()
- `abstol=1e-8`: the parameter of DifferentialEquations.ODEProblem.solve()

# Examples
KT simulation.

```julia-rep1
julia> K = 0.155
julia> T = 80.5
julia> target_δ_rad = 30.0 * π / 180.0
julia> target_Ψ_rad_deviation = 10.0 * π / 180.0
julia> time_second_interval = 0.01
julia> end_time_second = 500.0
julia> time_list = start_time_second:time_second_interval:end_time_second
julia> u_list, v_list, r_list, x_list, y_list, Ψ_list, δ_list, = kt_zigzag_test(
        K,
        T,
        time_list,
        target_δ_rad,
        target_Ψ_rad_deviation,
        u0 = 10.0, 
        x0 = 0.0,
        y0 = 0.0,
    )
```
"""
function kt_zigzag_test(
    K,
    T,
    time_list,
    target_δ_rad,
    target_Ψ_rad_deviation;
    u0=10.0,
    v0=0.0,
    r0=0.0,
    x0=0.0,
    y0=0.0,
    Ψ0=0.0,
    δ0=0.0,
    δ_rad_rate=10.0 * π / 180,
    algorithm=Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
)
    target_Ψ_rad_deviation = abs(target_Ψ_rad_deviation)

    final_δ_list = zeros(length(time_list))
    final_u_list = zeros(length(time_list))
    final_v_list = zeros(length(time_list))
    final_r_list = zeros(length(time_list))
    final_x_list = zeros(length(time_list))
    final_y_list = zeros(length(time_list))
    final_Ψ_list = zeros(length(time_list))

    next_stage_index = 1
    target_δ_rad = -target_δ_rad  # for changing in while loop
    Ψ = Ψ0
    while next_stage_index < length(time_list)
        target_δ_rad = -target_δ_rad
        start_index = next_stage_index

        # Make delta list
        δ_list = zeros(length(time_list) - start_index + 1)
        if start_index == 1
            δ_list[1] = δ0
            u0 = u0
            v0 = v0
            r0 = r0
            x0 = x0
            y0 = y0
        else
            δ_list[1] = final_δ_list[start_index-1]
            u0 = final_u_list[start_index-1]
            v0 = final_v_list[start_index-1]
            r0 = final_r_list[start_index-1]
            x0 = final_x_list[start_index-1]
            y0 = final_y_list[start_index-1]
        end

        for i in (start_index+1):length(time_list)
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

        u, v, r, x, y, Ψ_list, δ = kt_simulate(
            K,
            T,
            time_list[start_index:end],
            δ_list,
            u0=u0,
            v0=v0,
            r0=r0,
            x0=x0,
            y0=y0,
            Ψ0=Ψ,
            algorithm=algorithm,
            reltol=reltol,
            abstol=abstol,
        )

        # get finish index
        target_Ψ_rad = Ψ0 + target_Ψ_rad_deviation
        if target_δ_rad < 0
            target_Ψ_rad = Ψ0 - target_Ψ_rad_deviation
        end
        over_index = findfirst(e -> e > target_Ψ_rad, Ψ_list)
        if target_δ_rad < 0
            over_index = findfirst(e -> e < target_Ψ_rad, Ψ_list)
        end
        if isnothing(over_index)
            final_δ_list[start_index:end] = δ_list
            final_u_list[start_index:end] = u
            final_v_list[start_index:end] = v
            final_r_list[start_index:end] = r
            final_x_list[start_index:end] = x
            final_y_list[start_index:end] = y
            final_Ψ_list[start_index:end] = Ψ_list
            next_stage_index = length(time_list) # break
        else
            Ψ = Ψ_list[over_index]
            next_stage_index = over_index + start_index - 1
            final_δ_list[start_index:next_stage_index-1] = δ_list[begin:over_index-1]
            final_u_list[start_index:next_stage_index-1] = u[begin:over_index-1]
            final_v_list[start_index:next_stage_index-1] = v[begin:over_index-1]
            final_r_list[start_index:next_stage_index-1] = r[begin:over_index-1]
            final_x_list[start_index:next_stage_index-1] = x[begin:over_index-1]
            final_y_list[start_index:next_stage_index-1] = y[begin:over_index-1]
            final_Ψ_list[start_index:next_stage_index-1] = Ψ_list[begin:over_index-1]
        end
    end
    final_u_list,
    final_v_list,
    final_r_list,
    final_x_list,
    final_y_list,
    final_Ψ_list,
    final_δ_list
end

function estimate_kt_lsm(data::ShipData)
    time_vec = data.time
    r_vec = data.r
    δ_vec = data.δ

    spl_r = Spline1D(time_vec, r_vec)
    f = (spl, x) -> derivative(spl, x)
    dr_vec = f(spl_r, time_vec)
    estimate_kt_lsm(r_vec, dr_vec, δ_vec)
end

function estimate_kt_lsm(r_vec, dr_vec, δ_vec)
    A = hcat(δ_vec, r_vec)
    Θ = A \ dr_vec
    T = -1.0 / Θ[2]
    K = Θ[1] * T
    K, T
end

function estimate_kt_lsm_time_window_sampling(data::ShipData, window_size::Int)
    time_vec = data.time
    r = data.r
    spl_r = Spline1D(time_vec, r)
    f = (spl, x) -> derivative(spl, x)
    dr = f(spl_r, time_vec)
    δ = data.δ
    n_samples = length(time_vec) - window_size
    K_samples = zeros(n_samples)
    T_samples = zeros(n_samples)
    for i in 1:n_samples
        K_samples[i], T_samples[i] =
            estimate_kt_lsm(r[i:i+window_size], dr[i:i+window_size], δ[i:i+window_size])
    end
    K_samples, T_samples
end

function create_model_for_mcmc_sample_kt(
    data::ShipData;
    σ_r_prior_dist::Distribution=Chi(5),
    K_prior_dist::Distribution=Uniform(0.01, 10.0),
    T_prior_dist::Distribution=truncated(Normal(100.0, 50.0), 10.0, 200.0),
)
    time_obs = data.time
    r_obs = data.r
    δ_obs = data.δ

    # create sytem model
    spl_δ = Spline1D(time_obs, δ_obs)
    function KT!(dX, X, p, t)
        r, δ = X
        K, T = p
        dX[1] = dr = 1.0 / T * (-r + K * δ) # dr = 
        dX[2] = dδ = derivative(spl_δ, t) # dδ = 
    end
    u0 = [r_obs[1]; δ_obs[1]]
    K_start = 0.10
    T_start = 60.0
    p = (K_start, T_start)
    prob1 = ODEProblem(KT!, u0, (time_obs[1], time_obs[end]), p)

    # create probabilistic model
    Turing.@model function fitKT(time_obs, r_obs, prob1)
        σ_r ~ σ_r_prior_dist
        K ~ K_prior_dist
        T ~ T_prior_dist

        p = (K, T)
        prob = remake(prob1, p=p)
        sol = solve(prob, Tsit5())
        predicted = sol(time_obs)
        for i in 1:length(predicted)
            r_obs[i] ~ Normal(predicted[i][1], σ_r) # index number of r is 1
        end
    end

    return fitKT(time_obs, r_obs, prob1)
end
