# Demand function with symmetric, multiple cohorts
# p_i: the price of product of firm i (can be different from p[b])
# b: the cohort firm i belongs to
# p: NK-vector of prices across firms
# f: awareness set distribution
# params: model parameters
function demand_symmetric(p_i, b, p, f, params)
    @unpack cohorts, Γ_bar, Ω, q, σ, κ = params
    N = cohorts[1]
    b_bar = length(cohorts)

    f = reshape(f, Tuple(fill(0:N, b_bar)))
    demand_value = 0
    summand(n) = sum(Tuple(n) .* (p ./ q).^(-1/σ))
    for n in CartesianIndices(f)
        if (n[b] > 0)
            demand_value += n[b] * f[n] * (summand(n) + (p_i^(-1/σ) - p[b]^(-1/σ))/(q[b]^(-1/σ)))^(σ*(κ-1)-1)
        end
    end
    demand_value *= Γ_bar^(1-κ)*Ω*(q[b]^(1/σ))*(p_i^(-1/σ-1))
    return demand_value
end

# return type-irrelevant result as expected (for function signature consistency)
demand_symmetric(p_i, k, b, p, f, params) = demand_symmetric(p_i, b, p, f, params)

# Single summand for demand function with 2-type under a single cohort
# p_i: the price of product of firm i (can be different from p[k])
# k: the qualtity type firm i belongs to
# p: k_bar-vector of prices across firms
# f: awareness set distribution on 1:(N+1)
# params: model parameters
# n: awareness set in 1:(N+1)
# n_k: number of k-type firms awared such that n_k <= n
function demand_asymmetric_single_cohort_summand(p_i, k, p, f, params, n, n_k)
    @unpack cohorts, q, σ, κ, N_ks = params
    N = sum(N_ks)
    k_other = setdiff([1, 2], k)[1] # sends 1 to 2 and vice-versa
    n_1 = k == 1 ? n_k : N - n_k
    n_2 = N - n_1

    weight = (n_k / N_ks[k]) * (binomial(N_ks[k], n_k) * binomial(N_ks[k_other], n - n_k) / binomial(N, n))
    summands = [n_1; n_2] .* (p ./ q).^(-1/σ)
    return weight * (sum(summands) + (p_i^(-1/σ) - p[k]^(-1/σ)) / (q[k])^(-1/σ))^(σ*(κ-1)-1)
end

# Demand function with 2-type under a single cohort
# p_i: the price of product of firm i (can be different from p[k])
# k: the qualtity type firm i belongs to
# p: k_bar-vector of prices across firms
# f: awareness set distribution on 1:(N+1)
# params: model parameters
function demand_asymmetric_single_cohort(p_i, k, p, f, params)
    @unpack cohorts, Γ_bar, Ω, q, σ, κ, N_ks = params
    N = sum(N_ks)
    b_bar = length(cohorts)

    f = reshape(f, Tuple(fill(0:N, b_bar)))
    demand_value = 0

    # only first cohort
    for n in CartesianIndices(f)
        if (n[1] > 0)
            demand_value += f[n] * sum((n_k -> demand_asymmetric_single_cohort_summand(p_i, k, p, f, params, n[1], n_k)).(1:n[1]))
        end
    end
    return Γ_bar^(1-κ)*Ω*q[k]^(1/σ)*p_i^(-1/σ-1) * demand_value
end

# return cohort-irrelevant result as expected (for function signature consistency)
demand_asymmetric_single_cohort(p_i, k, b, p, f, params) = demand_asymmetric_single_cohort(p_i, k, p, f, params)

# Single summand for demand function with 2-type under multiple cohorts
# p_i: the price of product of firm i (can be different from p[k])
# k: the qualtity type firm i belongs to
# b: the cohort firm i belongs to
# p: k_bar by b_bar matrix of prices across qualities and cohorts
# f: awareness set distribution on 1:(N+1)
# params: model parameters
# n: k_bar-vector of whose each element is in 1:N
# n_k: number of k-type firms awared such that n_k <= n
function demand_asymmetric_summand(p_i, k, b, p, f, params, n, n_k)
    @unpack cohorts, q, σ, κ, N_ks = params
    b_bar = length(cohorts)
    N = sum(N_ks)

    k_other = setdiff([1, 2], k)[1] # sends 1 to 2 and vice-versa
    n_k1 = k == 1 ? Tuple(n_k) : N .- Tuple(n_k) # b_bar-vector representing 1st elements of n_k
    n_k2 = N .- n_k1 # b_bar-vector representing 2nd elements of n_k

    # given b_, compute the probability of choosing n_k[b_] (hypergeometric distribution)
    choose_n_kb_(b_) = (binomial(N_ks[k], n_k[b_]) * binomial(N_ks[k_other], n[b_] - n_k[b_]) /
                        binomial(N, n[b_]))
    # probability that firm i is awared and n_k firms are awared out of type k
    weight = (n_k[b] / N_ks[k]) * prod(choose_n_kb_.(1:b_bar))
    # sum on p/q over all cohorts and types
    summands = [n_k1 .* (p[1,:] ./ q[1]).^(-1/σ);
                n_k2 .* (p[2,:] ./ q[2]).^(-1/σ)]
    return weight * (sum(sum.(summands)) + (p_i^(-1/σ) - p[k,b]^(-1/σ)) / (q[k])^(-1/σ))^(σ*(κ-1)-1)
end


# Demand function with 2-type under multiple cohorts
# p_i: the price of product of firm i (can be different from p[k])
# k: the qualtity type firm i belongs to
# b: the cohort firm i belongs to
# p: k_bar by b_bar matrix of prices across qualities and cohorts
# f: awareness set distribution on 1:(N+1)
# params: model parameters
function demand_asymmetric(p_i, k, b, p, f, params)
    @unpack cohorts, Γ_bar, Ω, q, σ, κ, N_ks = params
    N = sum(N_ks)
    b_bar = length(cohorts)

    f = reshape(f, Tuple(fill(0:N, b_bar)))
    V_ns = CartesianIndices(Tuple(fill(0:N, b_bar))) # ⋃_n {V(n)}
    demand_value = 0
    for n in CartesianIndices(f)
        if (n[b] > 0)
            # extract V(n) by finding n_k in V_ns such that n_k <= n (elementwise) is true
            V_n = V_ns[(n_k -> n_k <= n).(V_ns)]
            # use v_k in V(n) such that v_k[b] > 0 (i.e., firm i in b cohort can be awared)
            V_n = V_n[(n_k -> n_k[b] > 0).(V_n)]
            demand_value += f[n] * sum((n_k -> demand_asymmetric_summand(p_i, k, b, p, f, params, n, n_k)).(V_n))
        end
    end
    return Γ_bar^(1-κ)*Ω*q[k]^(1/σ)*p_i^(-1/σ-1) * demand_value
end
