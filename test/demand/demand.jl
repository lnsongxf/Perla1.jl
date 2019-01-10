@testset "Demand functions" begin

    @testset "Monotonicity of demand functions" begin 
        for (params, demand) in ((params_demand_default_symmetric(), demand_symmetric),
            (params_demand_default_asymmetric_single_cohort(), demand_asymmetric_single_cohort),
            (params_demand_default_asymmetric(), demand_asymmetric))
            N = sum(params.N_ks)
            k_bar = length(params.N_ks)
            b_bar = length(params.cohorts)
            p = fill(1.0, k_bar, b_bar)
            f = fill(1/(N+1)^b_bar, (N+1)^b_bar)

            for k in 1:k_bar
                for b in 1:b_bar
                    # confirm that the demand function is decreasing in price
                    @test all(diff((p_i -> demand(p_i, k, b, p, f, params)).(1:100)) .<= 0)
                end
            end

        end
    end 
 
    @testset "multiple-cohort-two-type demand on single cohort nests = single-cohort-two-type demand" begin
        # asymmetric 1 cohort = asymmetric single cohort
        params = params_demand_default_asymmetric_single_cohort()
        N = sum(params.N_ks)
        b_bar = 1
        f = fill(1/(N+1), (N+1))

        b = 1 # fix cohort to one
        p = [1.0; 1.0] # fixed prices per-type
        for k in 1:2
            for p_i in 1.0:.1:3.0
                @test demand_asymmetric_single_cohort(p_i, k, p, f, params) ≈ 
                    demand_asymmetric(p_i, k, b, p, f, params) 
            end
        end

    end
end
