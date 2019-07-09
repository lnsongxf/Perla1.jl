# Compute Nash equilbria of prices; returns k_bar by b_Bar matrix of equilibrium prices
# f: awareness set distribution
# params: model parameters
# profit: profit function that takes (p_i, k, b, p, f, params)
function compute_price_equilibrium(f, params, profit; autodiff=:forward,
        iv = fill(params.mc, length(params.N_ks), length(params.cohorts)))
    ks = 1:length(params.N_ks)
    bs = 1:length(params.cohorts)
    profit_derivative = get_profit_derivative(profit)

    # return the gradient with respect to p (vectorized if p is a matrix)
    D_p = (p, f, params) -> vcat([profit_derivative(p[k,b], k, b, p, f, params) for k in ks, b in bs]...)

    sol = nlsolve(p -> D_p(p, f, params), iv, autodiff = autodiff)
    return (p = sol.zero, converged = converged(sol), solution = sol)
end
