# @testset "regression tests for industry_equilibrium" begin
#     params = stationary_params_default()
#     settings = stationary_settings_default()

#     # For a given industry
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, ), params)
#     θ_i_bounds_test(x, y, z) = [0.01; 0.2]

#     iv_test(params) = [0.07, 1.2]
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], μ_i_bounds = (x, y, z) -> [1.0, 1.0], iv = iv_test, θ_i_bounds = θ_i_bounds_test)
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.059999909404971456 atol = 1e-6
#     @test sol.f_converged

#     iv_test(params) = [0.05, 1.2]
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.05, 1.2)], μ_i_bounds = (x, y, z) -> [1.0, 1.0], iv = iv_test, θ_i_bounds = θ_i_bounds_test)
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.05
#     @test sol.zero[1] ≈ 0.059999914036147385 atol = 1e-6
#     @test sol.f_converged

#     iv_test(params) = [0.07, 1.3]
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], μ_i_bounds = (x, y, z) -> [1.0, 1.0], iv = iv_test, θ_i_bounds = θ_i_bounds_test)
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.05999991374646249 atol = 1e-6
#     @test sol.f_converged

#     # Completely different value
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, ν = 0.02), params)
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], μ_i_bounds = (x, y, z) -> [1.0, 1.0], iv = iv_test)
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.06415824644776288 atol = 1e-6
#     @test sol.f_converged
# end

# @testset "regression tests for industry_equilibrium under full setup" begin
#     iv_test(params) = [0.07, 1.2]
#     params = stationary_params_default()
#     settings = stationary_settings_default()

#     # η != 0.0
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, η = 0.02), params)
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], ignore_p_i_deviations = false, μ_i_bounds = (x, y, z) -> [1.0, 1.0], iv = iv_test)
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.06073385879817747 atol = 1e-6
#     @test sol.f_converged

#     # μ_i is not fixed
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, η = 0.1), params)
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], ignore_p_i_deviations = false, iv = iv_test, μ_i_bounds = ((x, y, z) -> [1.0; 1.1]))
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.06476724313165225 atol = 1e-6
#     @test sol.zero[2] ≈ 1.0220999258085721 atol = 1e-6
#     @test sol.f_converged

#     # μ_i is not fixed (lower cost for μ_i)
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, η = 0.05), params)
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], ignore_p_i_deviations = false, iv = iv_test, μ_i_bounds = ((x, y, z) -> [1.0; 1.1]))
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.06346214810200253 atol = 1e-6
#     @test sol.zero[2] ≈ 1.043617742891412 atol = 1e-6
#     @test sol.f_converged

#     # μ_i is not fixed (higher cost for μ_i)
#     industry_params = merge(industry_params_default(d = Perla1.d_alt, ϕ = 2.0, η = 0.6), params)
#     industry_settings = industry_settings_default(θ_i_μ_i_iv = (x, y, z) -> [(0.07, 1.2)], ignore_p_i_deviations = false, iv = iv_test, μ_i_bounds = ((x, y, z) -> [1.0; 1.1]))
#     sol = industry_equilibrium(industry_params, merge(settings, industry_settings))
#     @test sol.initial_x[1] ≈ 0.07
#     @test sol.zero[1] ≈ 0.10432590867444143 atol = 1e-6
#     @test sol.zero[2] ≈ 1.0047606937143156 atol = 1e-6
#     @test sol.f_converged
# end
