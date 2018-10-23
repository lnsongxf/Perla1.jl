@testset "solve_transition_dynamics" begin


    # simple duopoly
    N = 2
    T = 100.0
    μ = 0.0
    θ = 0.06 # baseline parameter from Perla16 (Appendix E.4)
    θ_d = 0.21 # baseline parameter from Perla16 (Appendix E.4)
    f0(a) = 0.5 # time-invariant f0
    f0_a(a) = (θ_d + θ) / (θ_d + θ * exp(θ_d + θ)*a) # (Appendix E.1)
    params_base = @with_kw (N = N, μ = μ, θ = θ, θ_d = θ_d, f0 = f0) # base parameters

    # time-invariant dynamics
    function Q_0!(df, f, p, t)
        @unpack N, μ, θ, θ_d, f0 = p
        df[1] = -(θ + θ_d*(1-f0(0.0)))*f[1] + μ*f[2]
        for i in 2:N
            df[i] = θ*((N+2-i)/N)*f[i-1] - (μ+θ*((N+1-i)/N))*f[i] + μ*f[i+1]
        end
        df[2] = (θ + θ_d*(1-f0(0.0)))*f[1] - (μ+θ*((N-1)/N))*f[2] + μ*f[3]
        df[end] = (θ/N)*f[N] - μ*f[N+1]
    end

    # time-variant dynamics
    function Q_a!(df, f, p, t)
        @unpack N, μ, θ, θ_d, f0 = p
        df[1] = -(θ + θ_d*(1-f0(t)))*f[1] + μ*f[2]
        for i in 2:N
            df[i] = θ*((N+2-i)/N)*f[i-1] - (μ+θ*((N+1-i)/N))*f[i] + μ*f[i+1]
        end
        df[2] = (θ + θ_d*(1-f0(t)))*f[1] - (μ+θ*((N-1)/N))*f[2] + μ*f[3]
        df[end] = (θ/N)*f[N] - μ*f[N+1]
    end

    # time values to be evaluated on for testing
    ts = range(0.0, stop = T, length = 50)

    # ===================================================================
    # functions for sanity check
    # ===================================================================
    # get Q in matrix form for sanity check
    function get_Q_matrix(N, μ, θ, θ_d, f0)
        dl = fill(μ, N)
        d = collect(-μ.-(N:-1:0).*(θ/N))
        du = collect((N:-1:1).*(θ/N))
        Q_basis = LinearAlgebra.Tridiagonal(dl,d,du)
        
        function Q(a)
            Q_basis[1,1] = -(θ + θ_d*(1-f0(a)))
            Q_basis[1,2] = θ + θ_d*(1-f0(a))
            return Q_basis
        end 

        return Q
    end

    # solve model with matrix Q
    function solve_transition_dynamics_matrix(Q_matrix, f_0, T)
        # solve transition dynamics given 
        # Q; N by N matrix generator
        # f_0; N vector of initial distribution
        # T; Float64 terminal time
        df(f,p,a) = Q_matrix(a)' * f
        prob = DifferentialEquations.ODEProblem(df,f_0,(0.0,T))
        return solve(prob);
    end

    # perform sanity check
    function sanity_check_dynamics(Q, params, f_0) 
        # solve dynamics
        sol_count = solve_transition_dynamics(Q, params, f_0, T)

        # solve dynamics, using matrix for benchmark
        Q_matrix = get_Q_matrix(params.N, params.μ, params.θ, params.θ_d, params.f0)
        sol_count_matrix = solve_transition_dynamics_matrix(Q_matrix, f_0, T)
        
        # average product awareness
        f_count(a) = dot(0:N, sol_count(a)) 
        f_count_matrix(a) = dot(0:N, sol_count_matrix(a)) # average awareness

        # check if f is a probability distribution for all t in (0, T)
        @test all(sum.(sol_count.(ts)) .≈ 1.0)

        # check if there is no forgetting, 
        # i.e., f_count is increasing.
        @test all(diff(f_count.(ts)) .> 0)

        # check if solutions are close to the ones based on matrix
        # (if Q is time-variant there are numerical errors, be more generous)
        f0s = (params.f0).(ts)
        is_Q_time_variant = any(y -> y != first(f0s), f0s)
        if (is_Q_time_variant)
            @test f_count.(0:0.1:T) ≈ f_count_matrix.(0:0.1:T) atol=1e+1
        else
            @test f_count.(0:0.1:T) ≈ f_count_matrix.(0:0.1:T) atol=1e-4
        end
    end

    # ===================================================================
    # unit test settings and executions
    # ===================================================================
    @testset "duopoly, time-invariant Q" begin
        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end

    @testset "duopoly, time-invariant Q, some firms recognized first" begin
        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [0.5; fill((1.0-0.5)/N, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end

    @testset "duopoly, time-variant Q" begin
        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a!, params, f_0)
    end

    @testset "duopoly, time-variant Q, some firms recognized first" begin
        # define generator
        params = params_base(N = N)
        # definte initial dist.
        f_0 = [0.5; fill((1.0-0.5)/N, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a!, params, f_0)
    end


    @testset "10 firms, time-invariant Q" begin
        N = 10 # 100 firms

        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end

    @testset "10 firms, time-invariant Q" begin
        N = 10 # 100 firms

        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end

    @testset "50 firms, time-invariant Q" begin
        N = 50 # 50 firms

        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end


    @testset "100 firms, time-invariant Q" begin
        N = 100 # 100 firms

        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end

    @testset "500 firms, time-invariant Q, some firms recognized first" begin
        N = 500

        # define generator
        params = params_base(N = N)

        # definte initial dist.
        f_0 = [0.5; fill((1.0-0.5)/N, N)]

        # perform sanity check
        sanity_check_dynamics(Q_0!, params, f_0)
    end
end
