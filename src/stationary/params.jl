function ℚf_default!(df, f, params, a)
    @unpack N, θ, θ_d, ζ = params

    # ℚ given a and f(a); baseline process from (C.1)
    df .= Tridiagonal(ζ*ones(N),
                [-(θ + θ_d * (1 - f[1])); -ζ.-((N-1):-1:0)/N*θ],
                [(θ + θ_d * (1 - f[1])); ((N-1):-1:1)/N*θ])' * f
end
# baseline settings and parameters for aggregate steady state
ss_iv_default(params) = [30.0; 60.0]
Ψ_default(μ, n, params) = μ + n - 1
stationary_params_default = @with_kw (N = 60, θ = 0.06, μ = 1.0, θ_d = 0.21,
    ζ = 0.0, σ = 0.15, κ = 3.5, ρ = 0.03, α = 0.28, δ_M = 0.056, δ_k = 0.07,
    z_M = 1.0, ℚf! = ℚf_default!, Ψ = Ψ_default)
stationary_settings_default = @with_kw (expectation_node_count = 20, ss_iv = ss_iv_default)

# baseline settings and parameters for industry equilibrium
θ_i_bounds_default(θ, μ, params) = [min(0.01,θ); max(0.2,θ)]
μ_i_bounds_default(θ, μ, params) = [min(μ, 1.0); max(3.0, μ + 0.1)]
p_i_bounds_default(p, mc, params) = ((1 + params.σ)*mc, p * 1.5)
iv_default(params) = [0.06, 1.5]

# CES d function (B.8)
function d_alt(θ_i, μ_i, params)
    @unpack η, ν, ϕ = params
    return 1/ν * ((1.0 - η )*θ_i^ϕ +  η *(μ_i-1)^ϕ  )^(2.0/ϕ)
end

function d_baseline(θ_i, μ_i, params)
    @unpack η, ν, ϕ = params
    return (θ_i^ϕ/(ν^(ϕ/2)) + abs(μ_i-1.0)^ϕ/η)^(2/ϕ)
end

# starting point for the optimizer (*NOT* fixedpoint) step
function θ_i_μ_i_iv_default(θ, μ, params)
  return [(θ, μ), (θ, μ * 1.3)]
end

industry_params_default = @with_kw (ν = 0.0178923, ϕ = 0.9, η = 0.0, d = d_baseline)
industry_settings_default = @with_kw (  θ_i_μ_i_iv = θ_i_μ_i_iv_default,
                                        μ_i_bounds = μ_i_bounds_default,
                                        θ_i_bounds = θ_i_bounds_default,
                                        p_i_bounds = p_i_bounds_default,
                                        ignore_p_i_deviations = true,
                                        clamp_θ_μ = true,
                                        iv = iv_default,
                                        ftol = 1E-6, g_tol = 1e-6,
                                        max_solver_iter = 10,
                                        max_iter = 1000,
                                        fixedpoint_beta = 1, fixedpoint_m = 5,
                                        trace_results = false,
                                        trace_iterations = false,
                                        Fminbox_algorithm = ConjugateGradient())
