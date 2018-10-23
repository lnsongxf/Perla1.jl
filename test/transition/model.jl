@testset "AwarenessModel" begin
    # simple duopoly 
    N = 2
    μ = 0.0
    θ = 0.06 # baseline parameter from Perla16 (Appendix E.4)
    θ_d = 0.0 # time-invariant transition matrix
    f0(a) = 0.5 # time-invariant transition matrix

    function LinearAlgebra.mul!(Y, Q::AwarenessModel, B::AbstractMatrix{T}, t) where {T}
        N = Q.N
        for j in 1:(N+1)
            y = Y[:,j]
            mul!(y, Q, B[:,j],t)
            Y[:,j] = y
        end
    end

    @testset "Check if duopoly model is valid" begin
        t = 0.0 # check on t = 0.0
        Q_expected = [-(θ + θ_d*(1-f0(t))) (θ + θ_d*(1-f0(t))) 0.0; 0.0 -θ/2 θ/2; 0.0 0.0 0.0]

        model_generated = AwarenessModel(N, μ, θ, θ_d, f0)

        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        expected = copy(eyes)
        generated = copy(eyes)
        LinearAlgebra.mul!(expected, Q_expected', eyes)
        LinearAlgebra.mul!(generated, model_generated, eyes, t)
        @test expected == generated
    end

    @testset "Check if the corresponding generator is a valid transition matrix" begin
        # Every row of the corresponding Q should sum up to 0;
        # since transpose of Q is used, need to check column
        t = 0.0 # check on t = 0.0

        # for N = 10
        N = 10
        model_generated = AwarenessModel(N, μ, θ, θ_d, f0)
        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        generated = copy(eyes)
        LinearAlgebra.mul!(generated, model_generated, eyes, t)
        for i in 1:(N+1)  # all column sums should be zero
            @test sum(generated[:,i]) ≈ 0.0 atol = 1e-8
        end

        # for N = 20, with non-zero θ_d
        N = 20
        θ_d = 0.15
        model_generated = AwarenessModel(N, μ, θ, θ_d, f0)
        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        generated = copy(eyes)
        LinearAlgebra.mul!(generated, model_generated, eyes, t)
        for i in 1:(N+1)  # all column sums should be zero
            @test sum(generated[:,i]) ≈ 0.0 atol = 1e-8
        end
    end
end
