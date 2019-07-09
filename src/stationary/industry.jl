# (B.4)
function Ψ_hat(p_i, p, μ_i, μ, n̂, params)
    Ψ(μ, n) = params.Ψ(μ, n, params)
    return Ψ(μ, n̂) + μ_i * (p_i / p)^(-1/params.σ) - μ
 end
# (B.7) value of the firm choosing μ_i when all other firms choose μ
function π_star_i(p_i, μ_i, μ, a, p, f_a, mc, Y, params)
    @unpack σ, κ, N = params
    E_n̂_m1 = E_by_n̂(n -> n*Ψ_hat(p_i, p, μ_i, μ, n, params)^(σ*(κ-1)-1), f_a)
    return μ_i*Y*(p - mc)*(1 - f_a[1])/N*p^(-κ) * E_n̂_m1
end
# (B.7) with uniform μ = μ_i = 1; value of each firm when all firms choose μ = 1
function π_star_i(μ, a, p, f_a, mc, Y, params)
    @unpack σ, κ, N = params
    Ψ(μ, n) = params.Ψ(μ, n, params)
    return Y*(p - mc)*(1 - f_a[1])/N*p^(-κ) * q(E_by_n̂(n -> Ψ(μ, n)^(σ*(κ-1)), f_a))  # (A.22) and # (A.7)
end
# (B.10) time 0 value of a firm
function v_0_i(θ_i, μ_i, θ, μ, a, p, f, mc, Y, quad_weights, params, settings)
   @unpack d, ν, η = params
   @unpack δ_M, ρ, N, θ_d, ζ = params
   r = δ_M + ρ

   if (μ_i ≈ μ) || settings.ignore_p_i_deviations
    p_i_val = p
   else
    p_i_val = p_i(p, μ_i, μ, a, f, mc, params, settings)
   end

   integrands = π_star_i.(p_i_val, μ_i, μ, a, p, f, mc, Y, Ref(params))
   return quad_weights ⋅ (θ_i / θ * integrands) - d(θ_i, μ_i, params) / N
end

# Equilibrium entrypoint
function industry_equilibrium(params, settings)
    ss_aggregate = steadystate(params, settings);
    r = params.δ_M + params.ρ;

    sol = fixedpoint(θ_μ -> industry_equilibrium_fixedpoint(θ_μ, ss_aggregate, params, settings), settings.iv(params),
                        beta = settings.fixedpoint_beta,
                        ftol = settings.ftol,
                        iterations = settings.max_iter)

    NLsolve.converged(sol) || throw(sol)
    # final clamp to recreate?
    return sol
end

# Equilibrium kernel (optimizer)
function industry_equilibrium_objective(θ_i_μ_i, params, settings, industry_aggregates)
    θ_i, μ_i = θ_i_μ_i
    @unpack θ, μ, a, p_a, f, mc, Y, pdv_weights = industry_aggregates
    return -v_0_i(θ_i, μ_i, θ, μ, a, p_a, f, mc, Y, pdv_weights, params, settings)
end

# Equilibrium kernel (fixedpoint)
function industry_equilibrium_fixedpoint(θ_μ, ss_aggregate, params, settings)
    θ, μ = θ_μ
    @unpack mc, Y, pdv_weights, a = ss_aggregate
    @unpack N, θ_d, ζ, δ_M, ρ = params
    θ_i_bounds = settings.θ_i_bounds(θ, μ, params);
    μ_i_bounds = settings.μ_i_bounds(θ, μ, params);
    r = params.δ_M + params.ρ

    industry_params = merge(params, (θ = θ, μ = μ,)) # swap out the industry ones
    @unpack E_n̂, E_n̂_m1, E_n̂_m2, f_0, f = solve_awareness_dynamics(industry_params, a; tstops = a)
    q_a = q.(E_n̂)
    p_a = Upsilon.(E_n̂_m1, E_n̂_m2, Ref(industry_params)) * mc
    industry_aggregates = (θ = θ, μ = μ, p_a = p_a, f = f, mc = mc, Y = Y, pdv_weights = pdv_weights, a = a)
    objective(θ_i_μ_i) = industry_equilibrium_objective(θ_i_μ_i, industry_params, settings, industry_aggregates)

    # switch depending on input
    if μ_i_bounds[1] ≈ μ_i_bounds[2]
        sol = Optim.optimize(x -> objective([x, μ]), θ_i_bounds[1], θ_i_bounds[2], Brent())
        Optim.converged(sol) || throw(sol)
        return vcat(sol.minimizer, [μ])
    else
        ivs = settings.θ_i_μ_i_iv(θ, μ, params)
        solutions = map(x -> Optim.optimize(objective, [θ_i_bounds[1], μ_i_bounds[1]], [θ_i_bounds[2], μ_i_bounds[2]], [x[1], x[2]],
                                Fminbox(settings.Fminbox_algorithm),
                                Optim.Options(g_tol = settings.g_tol,
                                            iterations = settings.max_solver_iter,
                                            show_trace = settings.trace_results)), ivs)
        # extract best value and return
        filter!(Optim.converged, solutions) # drop any non-converged solutions
        length(solutions) > 0 || error("No solutions converged.") # if all failed to converge, then error out
        minima = [sol.minimum for sol in solutions];
        minimizers = [sol.minimizer for sol in solutions];
        return minimizers[findmin(minima)[2]] # return the minimizer leading to the best point
    end
end

# =================================
# off equilibrium pricing equations
# =================================
function p_i_root(p_i, p_a, μ_i, μ, f_a, mc, params)
    # implementation of (B.3)
    # calls g_i and D_g_i
    @unpack σ = params
    (-p_i + mc + σ*mc)*g_i(p_i, p_a, μ_i, μ, f_a, mc, params)+p_i*(p_i-mc)σ*D_g_i(p_i, p_a, μ_i, μ, f_a, mc, params)
end

function g_i(p_i, p_a, μ_i, μ, f_a, mc, params)
    # implementation of (B.5)  where p_a scalar
    # f_a is a vector here, at the particular a.
    @unpack σ, κ = params
    E_n̂_m1 = E_by_n̂(n -> n*Ψ_hat(p_i, p_a, μ_i, μ, n, params)^(σ*(κ-1)-1), f_a)
    return p_a^(1/σ-κ+1) * E_n̂_m1
end

function D_g_i(p_i, p_a, μ_i, μ, f_a, mc, params)
    # implementation of (B.6)  where p_a scalar
    # f_a is a vector here, at the particular a.
    @unpack σ, κ = params
    E_n̂_m2 = E_by_n̂(n -> n*Ψ_hat(p_i, p_a, μ_i, μ, n, params)^(σ*(κ-1)-2), f_a)
    return μ_i*(1-σ*(κ-1)) / σ * p_a^(1/σ-κ) * (p_i / p_a)^(-1-1/σ) * E_n̂_m2
end

function p_i(p, μ_i, μ, a, f, mc, params, settings)
    # compute the price from (B.3) for every a and return a vector
    # p is a vector of prices for every a
    # f is the evolution from the equilibrium.  Use it for calculating expectations
    p_i = similar(p)
    for a_index in eachindex(a)
        f_a = f[a_index]
        p_a = p[a_index]
        p_i[a_index] = find_zero(p_i -> p_i_root(p_i, p_a, μ_i, μ, f_a, mc, params),
                        settings.p_i_bounds(p_a, mc, params))
    end

    return p_i
end

# ======================================================
# auxiliary functions, testing/inspection functions, etc.
# ======================================================

function industry_aggregates_for_tests(θ, μ, params, settings)
    ss_aggregate = steadystate(params, settings);
    @unpack mc, Y, pdv_weights, a = ss_aggregate
    @unpack N, θ_d, ζ, δ_M, ρ = params
    industry_params = merge(params, (θ = θ, μ = μ,)) # swap out the industry ones
    @unpack E_n̂, E_n̂_m1, E_n̂_m2, f_0, f = solve_awareness_dynamics(industry_params, a; tstops = a)
    q_a = q.(E_n̂)
    p_a = Upsilon.(E_n̂_m1, E_n̂_m2, Ref(industry_params)) * mc
    industry_aggregates = (θ = θ, μ = μ, p_a = p_a, f = f, mc = mc, Y = Y, pdv_weights = pdv_weights, a = a)
    return (industry_aggregates = industry_aggregates, ss_aggregate  = ss_aggregate) # can feed these directly to industry_equilibrium_objective
end
