@testset "AwarenessModel" begin
    # simple duopoly 
    N = 2
    μ = 0.0
    θ = 0.06 # baseline parameter from Perla16 (Appendix E.4)
    θ_d = 0.0 # time-invariant transition matrix
    f0(a) = 0.5 # time-invariant transition matrix

    cohorts = (N, )
    # base parameters
    params_base = @with_kw (μ = μ, θ = θ, θ_d = θ_d, f0 = f0, cohorts = cohorts) 

    # time-variant dynamics, with single cohort
    function Q_a!(df, f, p, t)
        @unpack μ, θ, θ_d, f0, cohorts = p
        N_1 = cohorts[1]
        e_1 = CartesianIndex((1,))
        current_cohort = 1 # single cohort case

        f = reshape(f, (cohorts.+1))
        df = reshape(df, (cohorts.+1))

        df[1] = -(θ + θ_d*(1-f0(t)))*f[1] + μ*f[2]
        for i in CartesianIndices(f)
            if (i[current_cohort] > 1 && i[current_cohort] <= N_1)
                i_previous = i - e_1
                i_forward = i + e_1
                df[i] = θ*((N+2-i[current_cohort])/N)*f[i_previous] - 
                        (μ+θ*((N+1-i[current_cohort])/N))*f[i] + 
                        μ*f[i_forward]
            end
        end
        df[2] = (θ + θ_d*(1-f0(t)))*f[1] - (μ+θ*((N-1)/N))*f[2] + μ*f[3]
        df[end] = (θ/N)*f[N] - μ*f[N+1]

        f = reshape(f, (N_1+1,))   
        df = reshape(df, (N_1+1,))
    end

    # helper function to apply dynamics per column
    function dynamics_by_col(Y, Q!, B::AbstractMatrix{T}, p, t) where {T}
        N = p.cohorts[1]
        for j in 1:(N+1)
            y = Y[:,j]
            Q!(y, B[:,j], p, t)
            Y[:,j] = y
        end
    end

    @testset "Check if duopoly model is valid" begin
        params = params_base()
        t = 0.0 # check on t = 0.0
        Q_expected = [-(θ + θ_d*(1-f0(t))) (θ + θ_d*(1-f0(t))) 0.0; 0.0 -θ/2 θ/2; 0.0 0.0 0.0]

        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        expected = copy(eyes)
        generated = copy(eyes)
        LinearAlgebra.mul!(expected, Q_expected', eyes)
        dynamics_by_col(generated, Q_a!, eyes, params, t)
        @test expected == generated
    end

    @testset "Check if the corresponding generator is a valid transition matrix" begin
        # Every row of the corresponding Q should sum up to 0;
        # since transpose of Q is used, need to check column
        t = 0.0 # check on t = 0.0

        # for N = 10
        N = 10
        params = params_base(cohorts = (N, ))
        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        generated = copy(eyes)
        dynamics_by_col(generated, Q_a!, eyes, params, t)
        for i in 1:(N+1)  # all column sums should be zero
            @test sum(generated[:,i]) ≈ 0.0 atol = 1e-8
        end

        # for N = 20, with non-zero θ_d
        N = 20
        θ_d = 0.15
        params = params_base(cohorts = (N, ), θ_d = θ_d)
        eyes = LinearAlgebra.Matrix{Float64}(I, N+1, N+1)
        generated = copy(eyes)
        dynamics_by_col(generated, Q_a!, eyes, params, t)
        for i in 1:(N+1)  # all column sums should be zero
            @test sum(generated[:,i]) ≈ 0.0 atol = 1e-8
        end
    end
end
