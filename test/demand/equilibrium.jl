@testset "Equilibrium solvers" begin
    @testset "Equilibrium prices are above mc" begin 
        for (params, demand) in ((params_demand_default_symmetric(), demand_symmetric),
            (params_demand_default_asymmetric_single_cohort(), demand_asymmetric_single_cohort),
            (params_demand_default_asymmetric(), demand_asymmetric))
            N = sum(params.N_ks)
            k_bar = length(params.N_ks)
            b_bar = length(params.cohorts)
            p = fill(1.0, k_bar, b_bar)
            f = fill(1/(N+1)^b_bar, (N+1)^b_bar)

            prices_above_mc = range(params.mc, stop = params.mc + 5, length=20)
            prices_below_mc = range(params.mc/2, stop=params.mc, length=20)

            for k in 1:k_bar
                for b in 1:b_bar
                    profit = get_profit(demand)
                    equilibrium_prices = compute_price_equilibrium(f, params, profit).p
                    @test all(equilibrium_prices .>= params.mc)
                end
            end
        end
    end
end
