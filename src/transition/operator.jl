e(i, N) = CartesianIndex((1:N .== i)*1...) # N-dimensional unit tuple at position i

# transition operator
function transition_operator_base!(df, f, p, t)
    # unpack params
    @unpack θ, cohorts, ts_cohort, ζ = p
    @assert length(unique(cohorts)) == 1 # assume that N_t is invariant across all t
    N = cohorts[1]
    K = length(cohorts) # number of cohorts
    setsize = (N+1)^K

    f = reshape(f, Tuple(fill(0:N, K)))
    df = reshape(df, Tuple(fill(0:N, K)))
    θ_t = θ(t)
    current_cohort = sum(t .> ts_cohort) + 1
    N_total = current_cohort*N;

    for i in CartesianIndices(f)
        # (outflow)
        full_aware_ratio = sum(Tuple(i) .== N) / current_cohort # ratio of cohorts fully awared
        any_aware_count = sum(Tuple(i) .!= 0) # number of cohorts that aren't zero
        df[i] = -((((N_total - sum(i)) / N_total) * θ_t) *
            (1-full_aware_ratio) + # (1-full_aware_ratio) cohorts aren't full (can be destinations)
            ζ) * f[i]

        # (inflow)
        for cohort in 1:current_cohort
            # if the current status has some awareness in a cohort
            # then there's in-flow from past in the cohort
            if (i[cohort] > 0)
                i_past = i - e(cohort, K)
                df[i] += ((N_total + 1 - sum(i)) / N_total) * (θ_t / current_cohort) * f[i_past]
            end
            # if the current status has no full-awareness in a cohort
            # then there's in-flow from more-awared status from the cohort source by forgetting
            if (i[cohort] < N)
                # Inflow amount depends on the current status:
                # 1.1.2 can be from 2.1.2 and 1.2.2;
                # note that the inflow from each case is ζ/3
                # 1.1.1 can be from 2.1.1 and 1.2.1 and 1.1.2
                # note that the inflow from each case is ζ/3
                # 0.0.1 can be from 1.0.1 and 0.1.1 and 0.0.2
                # note that the inflow from case 1 is ζ/2 while case 2/3 is ζ
                # from source state forgetting is dispersed by any_aware_count_in_source_state
                # which is (any_aware_count + (i[cohort] == 0)) in the current state because
                # if i[cohort] == 0 then the source state has any_aware_count + 1 (1 from cohort)
                # if i[cohort] != 0 then the source state has any_aware_count
                # (since the source state has i[cohort] > 0)
                i_future = i + e(cohort, K)
                df[i] += ζ / (any_aware_count + (i[cohort] == 0)) * f[i_future]
            end
        end
    end
    df[1] += ζ * f[1] # there is no forgetting in the state where no product is recognized
end

# auxiliary constructors
get_transition_operator(params) = get_transition_operator(transition_operator_base!, params)
get_transition_operator(transition_operator_base!, params) = MatrixFreeOperator(transition_operator_base!, (params, 0.0),
    size = (set_size(params), set_size(params)), opnorm=p->0.1)
