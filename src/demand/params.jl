# scalars used in defaults for all cases
Ω = 1.0
σ = 0.15
mc = 1.0
κ = 3.5
Γ_bar(σ, κ) = gamma(1-σ*(κ-1))^(1/(1-κ))
q_min = 2.0
q_max = 3.0

params_demand_default_base = @with_kw (Ω = Ω, σ = σ, mc = mc, κ = κ, Γ_bar = Γ_bar(σ, κ))
# multiple cohorts on symmetric firms
params_demand_default_symmetric = @with_kw (Ω = Ω, σ = σ, mc = mc, κ = κ, Γ_bar = Γ_bar(σ, κ), N = 4, N_ks = (N,),
                                            cohorts = (N, N), q = range(q_min, q_max, length = length(cohorts)))
# single cohort on binary firm types
params_demand_default_asymmetric_single_cohort = @with_kw (Ω = Ω, σ = σ, mc = mc, κ = κ, Γ_bar = Γ_bar(σ, κ),
                                                            N_k = 4, N = 2*N_k, N_ks = (N_k, N_k), cohorts = (N,),
                                                            q = range(q_min, q_max, length = length(N_ks)))
# multiple cohorts on binary firm types
params_demand_default_asymmetric = @with_kw (Ω = Ω, σ = σ, mc = mc, κ = κ, Γ_bar = Γ_bar(σ, κ),
                                            N_k = 4, N = 2*N_k, N_ks = (N_k, N_k), cohorts = (N, N),
                                            q = range(q_min, q_max, length = length(N_ks)))
