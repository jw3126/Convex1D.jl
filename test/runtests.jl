using Convex1D
using Test

@testset "minimize" begin
    res = @inferred minimize(x -> x^2, (-1.1, 1), atol=1e-3)
    @test res.minimizer ≈ 0 atol=1e-3
    @test res.minimum === res.minimizer^2

    for _ in 1:100
        x_opt = randn()
        f_opt = randn()
        s = 10rand()
        f = x -> s*(x - x_opt)^2 + f_opt
        atol = eps(x_opt)^(1/3)
        res = minimize(f, (x_opt - 10rand(), x_opt + 10rand()), atol=atol/10)

        @test x_opt ≈ res.minimizer atol=atol
        @test f_opt ≈ res.minimum
        @test res.nevals < 100 # quick convergence

        x_res, f_res = minimize(f, (x_opt, x_opt + 10rand()))
        @test x_opt ≈ x_res atol=atol
        @test f_res ≈ f_res atol=atol

        x_res, f_res = minimize(f, (x_opt-10rand(), x_opt), atol=atol)
        @test x_opt ≈ x_res atol=atol
        @test f_res ≈ f_res
    end
end
