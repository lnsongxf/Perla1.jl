ts_cohort_default_two_cohorts = [5.0]
θ_default_two_cohorts(t) = t < ts_cohort_default_two_cohorts[end] ? 0.4 : 0.6
params_transition_default_two_cohorts = @with_kw (θ = θ_default_two_cohorts,
    cohorts = (2, 2), ts_cohort = ts_cohort_default_two_cohorts, μ = 0.1)

ts_cohort_default_three_cohorts = [5.0; 10.0]
θ_default_three_cohorts(t) = t < ts_cohort_default_three_cohorts[1] ? 0.4 : (t < ts_cohort_default_three_cohorts[end] ? 0.6 : 1.2)
params_transition_default_three_cohorts = @with_kw (θ = θ_default_three_cohorts,
    cohorts = (2, 2, 2), ts_cohort = ts_cohort_default_three_cohorts, μ = 0.1)

ts_cohort_default_four_cohorts = [5.0; 10.0; 15.0]
θ_default_four_cohorts(t) = t < ts_cohort_default_four_cohorts[1] ? 0.4 : (t < ts_cohort_default_four_cohorts[2] ? 0.6 : (t < ts_cohort_default_four_cohorts[3] ? 1.2 : 1.5))
params_transition_default_four_cohorts = @with_kw (θ = θ_default_four_cohorts,
    cohorts = (2, 2, 2, 2), ts_cohort = ts_cohort_default_three_cohorts, μ = 0.1)

set_size(params) = prod(params.cohorts.+ 1)
