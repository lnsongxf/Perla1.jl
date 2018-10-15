@testset "solve_transition_dynamics" begin
    # simple duopoly
    N = 2
    T = 100.0
    μ = 0.0
    θ = 0.06 # baseline parameter from Perla16 (Appendix E.4)
    θ_d = 0.21 # baseline parameter from Perla16 (Appendix E.4)

    # time values to be evaluated on for testing
    ts = range(0.0, stop = T, length = 50)

    # function for sanity check
    function sanity_check_dynamics(Q, f_0) 
        # solve dynamics
        sol_count = solve_transition_dynamics(Q, f_0, T)

        # average product awareness
        f_count(a) = dot(0:N, sol_count(a)) 

        # check if f is a probability distribution for all t in (0, T)
        @test all(sum.(sol_count.(ts)) .≈ 1.0)

        # check if there is no forgetting, 
        # i.e., f_count is increasing.
        @test all(diff(f_count.(ts)) .> 0)
    end

    @testset "duopoly, time-invariant Q" begin
        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)
        Q(a) = Q_a(0) # time invariant dynamics

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q, f_0)
    end

    @testset "duopoly, time-variant Q" begin
        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end

    @testset "duopoly, time-variant Q, some firms recognized first" begin
        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [0.5; fill((1.0-0.5)/N, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end

    @testset "100 firms, time-variant Q" begin
        N = 100 # 100 firms

        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end

    @testset "100 firms, time-variant Q" begin
        N = 100 # 100 firms

        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end

    @testset "10000 firms, time-variant Q" begin
        N = 10000

        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [1.0; fill(0.0, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end

    @testset "10000 firms, time-variant Q, some firms recognized first" begin
        N = 10000

        # define generator
        Q_a = get_Q(N, μ, θ, θ_d)

        # definte initial dist.
        f_0 = [0.5; fill((1.0-0.5)/N, N)]

        # perform sanity check
        sanity_check_dynamics(Q_a, f_0)
    end


    
end
