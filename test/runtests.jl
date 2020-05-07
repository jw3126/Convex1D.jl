using Convex1D
using Test

@testset "minimize" begin
    res = @inferred minimize(x -> x^2, (-1.1, 1), atol=1e-3)
    @test res.minimizer ≈ 0 atol=1e-3
    @test res.minimum === res.minimizer^2

    x = [1,2]
    y = [2,4]
    f = c -> sum(abs, c*x - y)
    sol = minimize(f, (-10, 10), atol=1e-10)
    @test sol.minimizer ≈ 2 atol=1e-10

    for _ in 1:100
        x_opt = randn()
        f_opt = randn()
        s = 0.1 + 10rand()
        p = 1 + rand()
        f = x -> s*abs(x - x_opt)^p + f_opt
        @assert f(x_opt) ≈ f_opt
        atol = eps(x_opt)^(1/3)
        res = minimize(f, (x_opt - 10rand(), x_opt + 10rand()), atol=atol)
        @test res.nevals < 100 # quick convergence
        @test res.xbounds[1] <= res.minimizer <= res.xbounds[2]
        @test res.xbounds[2] - res.xbounds[1] <= atol
        # @test res.ybounds[1] <= res.minimum <= res.ybounds[2]
        # @test res.ybounds[1] <= f_opt <= res.ybounds[2]

        x_res, f_res = res
        @test x_opt ≈ x_res atol=atol
        @test f(x_res) == f_res

        res = minimize(f, (x_opt, x_opt + 10rand()), atol=atol)
        x_res, f_res = res
        @test x_opt ≈ x_res atol=atol
        @test f(x_res) == f_res

        I = (x_opt-10rand(), x_opt)
        res = minimize(f, I, atol=atol)
        x_res, f_res = res
        if !(isapprox(x_opt, x_res, atol=atol))
            @show x_opt
            @show f(x_opt)
            @show f_opt
            @show s
            @show p
            @show res
            @show I
            @show atol
            error()
        end
        @test x_opt ≈ x_res atol=atol
        @test f(x_res) == f_res
    end
end
