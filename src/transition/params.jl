ts_cohort_default_two_cohorts = [5.0]
θ_default_two_cohorts(t) = t < ts_cohort_default_two_cohorts[end] ? 0.4 : 0.6
params_transition_default_two_cohorts = @with_kw (θ = θ_default_two_cohorts, 
    cohorts = (2, 2), ts_cohort = ts_cohort_default_two_cohorts, μ = 0.1)

# ts_cohort_default_three_cohorts = [5.0; 10.0]