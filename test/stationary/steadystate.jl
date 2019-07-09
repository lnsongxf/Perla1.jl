@testset "regression tests for steadystate_exogenous μ = 1" begin
    params = stationary_params_default()
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # regresssion tests
    @test ss.Q ≈ 0.7546898959261538
    @test ss.B ≈ 0.802251370104198
    @test ss.k ≈ 38.898217160871326
    @test ss.M ≈ 64.61497867287925
    @test ss.C ≈ 7.550906407714338
    @test ss.capital_share ≈ 0.22463038362917545
    @test ss.labor_share ≈ 0.5776209864750226
    @test ss.profit_share ≈ 0.197748629895802
    @test ss.mc ≈ 3.998543431923668
    @test ss.Z ≈ 4.984152824075987
    @test ss.w ≈ 8.024438060242373
    @test ss.Y ≈ 13.892220414656569
    @test ss.TobinsQ ≈ 1.077269705336297
    @test ss.p[1] ≈ 5.5391113031169334
    @test ss.p[end] ≈ 4.60745274893
    @test ss.a[1] ≈ 1.2596408873569376
    @test ss.a[end] ≈ 1187.9360093859957
    @test ss.q[1] ≈ 1.0106460062423956
    @test ss.q[end] ≈ 4.054153800046078

    # check if elements of f at each a sum up to one
    for f_a in ss.f 
        @test sum(f_a) ≈ 1.0
    end
end

@testset "regression tests for steadystate_exogenous, when μ > 1" begin
    params = stationary_params_default(μ = 1.1)
    settings = stationary_settings_default()
    ss = steadystate(params, settings)

    # regresssion tests
    @test ss.Q ≈ 0.7694677464576509
    @test ss.B ≈ 0.8098122639789883
    @test ss.k ≈ 40.13443977752665
    @test ss.M ≈ 66.66850461503046
    @test ss.C ≈ 7.790881449318419
    @test ss.capital_share ≈ 0.22674743391411675
    @test ss.labor_share ≈ 0.5830648300648715
    @test ss.profit_share ≈ 0.1901877360210117
    @test ss.mc ≈ 4.1281807649266
    @test ss.Z ≈ 5.0977009716337305
    @test ss.w ≈ 8.357492967493016
    @test ss.Y ≈ 14.33372849218699
    @test ss.TobinsQ ≈ 1.074315307918037
    @test ss.p[1] ≈ 5.521011376538905
    @test ss.p[end] ≈ 4.75680844967857
    @test ss.a[1] ≈ 1.2596408873569376
    @test ss.a[end] ≈ 1187.9360093859957
    @test ss.q[1] ≈ 1.0465906204405908
    @test ss.q[end] ≈ 4.0577974567349

    # check if elements of f at each a sum up to one
    for f_a in ss.f 
        @test sum(f_a) ≈ 1.0
    end
end