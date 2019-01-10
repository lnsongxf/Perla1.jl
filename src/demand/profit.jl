# given demand function, return a profit function that takes (p_i, k, b, p, f, params)
get_profit(demand) = (p_i, k, b, p, f, params) -> ((p_i - params.mc) * demand(p_i, k, b, p, f, params))
