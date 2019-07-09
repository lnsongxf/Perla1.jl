@testset "p_i, when μ_i = μ = 1" begin
    # solve steady states
    params = stationary_params_default()
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # extract p, a, f, mc
    @unpack p, a, f, mc = ss 

    # identical μ
    μ = params.μ
    μ_i = μ

    # expect identical pricing
    settings = industry_settings_default()
    p_i_off = p_i(p, μ_i, μ, a, f, mc, params, settings)
    @test p ≈ p_i_off
end

@testset "p_i, when μ_i > μ = 1" begin
    # solve steady states
    params = stationary_params_default()
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # extract p, a, f, mc
    @unpack p, a, f, mc = ss 

    # higher μ_i
    μ = params.μ
    μ_i = μ + 0.1

    # expect higher pricing
    settings = industry_settings_default()
    p_i_off = p_i(p, μ_i, μ, a, f, mc, params, settings)
    @test all(p .< p_i_off)

    # regression tests
    @test p_i_off[1] ≈ 5.540222769066783
    @test p_i_off[3] ≈ 5.091922516464228
    @test p_i_off[end] ≈ 4.608343292733195
end