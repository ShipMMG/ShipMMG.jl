function kt_simulate(
    time_list,
    δ_list,
    K,
    T,
    u0,
    x0 = 0.0,
    y0 = 0.0,
    ψ0 = 0.0,
    r0 = 0.0,
    algorithm = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-8,
)

    spl_δ = Spline1D(time_list, δ_list)

    function kt_eom!(dX, X, KT, t)
        u, r, δ = X
        K, T = KT
        dX[1] = du = 0.0
        dX[2] = dr = 1.0 / T * (-r + K * δ)
        dX[3] = dδ = derivative(spl_δ, t)
    end

    X0 = [u0; r0; δ_list[1]]
    KT = [K, T]
    prob = ODEProblem(kt_eom!, X0, (time_list[1], time_list[end]), KT)
    sol = solve(prob, algorithm, reltol = reltol, abstol = abstol)

    results = hcat(sol.u...)
    time = sol.t
    u = results[1, :]
    r = results[2, :]
    δ = results[3, :]
    time, u, r, δ
end
