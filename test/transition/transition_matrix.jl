@testset "get_Q" begin
    # simple duopoly 
    N = 2
    μ = 0.0
    θ = 0.5
    θ_d = 0.00

    f_0(a) = (θ_d + θ) / (θ_d + θ*exp((θ_d + θ)*a))
    Q_expected(a) = [-(θ + θ_d*(1-f_0(a))) (θ + θ_d*(1-f_0(a))) 0.0; 0.0 -θ/2 θ/2; 0.0 0.0 0.0]

    Q_generated = get_Q(N, μ, θ, θ_d)

    @test Q_generated.(-3:0.01:3) == Q_expected.(-3:0.01:3)
end
