@testset "π_star_i, when μ_i = μ = 1, p_i = p" begin
    # solve steady states
    params = stationary_params_default()
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # extract p, a, f, f_0, mc, Y
    @unpack p, a, f, f_0, mc, Y = ss 

    # higher μ_i
    μ = params.μ
    μ_i = μ
    p_i_off = p

    # expect identical profit
    π_star = π_star_i.(μ, a, p, f, mc, Y, Ref(params))
    π_star_off = π_star_i.(p_i_off, μ_i, μ, a, p, f, mc, Y, Ref(params))
    @test π_star ≈ π_star_off
end

@testset "π_star_i, when μ_i > μ = 1, p_i = p" begin
    # solve steady states
    params = stationary_params_default()
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # extract p, a, f, f_0, mc, Y
    @unpack p, a, f, f_0, mc, Y = ss 

    # higher μ_i
    μ = params.μ
    μ_i = μ + 0.1
    p_i_off = p

    # expect higher profit
    π_star = π_star_i.(μ, a, p, f, mc, Y, Ref(params))
    π_star_off = π_star_i.(p_i_off, μ_i, μ, a, p, f, mc, Y, Ref(params))
    @test all(π_star .< π_star_off)

    # regression tests
    @test p_i_off[1] ≈ 5.5391113031169334
    @test p_i_off[5] ≈ 4.759693946483811
    @test p_i_off[end] ≈ 4.60745274893
end