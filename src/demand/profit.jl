# given demand function, return a profit function that takes (p_i, k, b, p, f, params)
get_profit(demand) = (p_i, k, b, p, f, params) -> ((p_i - params.mc) * demand(p_i, k, b, p, f, params))
# given profit function, return a profit function that takes (p_i, k, b, p, f, params)
get_profit_derivative(profit) = (p_i, k, b, p, f, params) -> derivative(p_i -> profit(p_i, k, b, p, f, params), p_i)