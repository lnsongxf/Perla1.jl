# default parameters for all cohort/quality cases
params_demand_default_base = @with_kw (Ω = 1.0, σ = 0.15, κ = 3.5, mc = 1.0, Γ_bar = gamma(1-σ*(κ-1))^(1/(1-κ))) 
# default parameters for multiple cohorts on symmetric firms
params_demand_default_symmetric = @with_kw (Ω = 1.0, σ = 0.15, κ = 3.5, mc = 1.0, Γ_bar = gamma(1-σ*(κ-1))^(1/(1-κ)),
    N = 4, N_ks = (N, ), cohorts = (N, N, ),
    q = range(2.0, 3.0, length = length(cohorts))) # this represents quality per cohort (symmetric within cohort)
# default parameters for a single cohort on two type firms
params_demand_default_asymmetric_single_cohort = @with_kw (Ω = 1.0, σ = 0.15, κ = 3.5, mc = 1.0, Γ_bar = gamma(1-σ*(κ-1))^(1/(1-κ)),
    N_k = 4, N = 2*N_k, N_ks = (N_k, N_k), cohorts = (N, ),
    q = range(2.0, 3.0, length = length(N_ks))) # this represents quality per type (asymmetric within cohort, symmetric within type)
# default parameters for multiple cohorts on two type firms
params_demand_default_asymmetric = @with_kw (Ω = 1.0, σ = 0.15, κ = 3.5, mc = 1.0, Γ_bar = gamma(1-σ*(κ-1))^(1/(1-κ)),
    N_k = 4, N = 2*N_k, N_ks = (N_k, N_k), cohorts = (N, N, ),
    q = range(2.0, 3.0, length = length(N_ks)) ) # this represents quality per type (asymmetric within cohort, symmetric within type)
