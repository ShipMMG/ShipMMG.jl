@testset "calc_position func" begin
    time_vec = [0.0, 2.0]
    u_vec = [1.0, 1.0]
    v_vec = [0.0, 0.0]
    r_vec = [0.0, 0.0]
    x, y, ψ = calc_position(time_vec, u_vec, v_vec, r_vec, x0=1.0, y0=1.0, ψ0=0.0)
    @test x[2] ≈ 3.0
    @test y[2] ≈ 1.0
    @test ψ[2] ≈ 0.0
end
