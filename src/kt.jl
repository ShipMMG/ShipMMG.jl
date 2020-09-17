function kt_simulate(time_list, δ_list, K, T, u0, x0=0.0, y0=0.0, ψ0=0.0, r0=0.0, algorithm = Tsit5(), reltol=1e-8,abstol=1e-8)

    spl_δ = Spline1D(time_list, δ_list)

    function kt_eom!(dX, X, KT, t)
        x, y, ψ, u, r, δ = X
        K, T = KT
        dX[1] = dx = u * cos.(ψ)
        dX[2] = dy = u * sin.(ψ)
        dX[3] = dψ = r
        dX[4] = du = 0.0
        dX[5] = dr = 1.0 / T * ( -r + K * δ)
        dX[6] = dδ = derivative(spl_δ, t)
    end

    X0 = [x0; y0; ψ0; u0; r0; δ_list[1]]
    KT = [K, T]
    prob = ODEProblem(kt_eom!, X0, (time_list[1], time_list[end]), KT)
    sol = solve(prob, algorithm, reltol=reltol,abstol=abstol)
    sol
end