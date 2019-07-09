module Perla1

using SparseArrays, DifferentialEquations, LinearAlgebra, ExponentialUtilities, DiffEqOperators, Parameters, SpecialFunctions, OffsetArrays, NLsolve
using Distributions, Expectations, Optim, Roots
using ForwardDiff: jacobian, derivative
import Parameters: @with_kw, @unpack

include("transition/solve-awareness-dynamics.jl")
include("transition/params.jl")
include("transition/operator.jl")
include("stationary/params.jl")
include("stationary/steadystate.jl")
include("stationary/industry.jl")
include("demand/params.jl")
include("demand/demand.jl")
include("demand/profit.jl")
include("demand/equilibrium.jl")

export solve_transition_dynamics, solve_awareness_dynamics
export params_demand_default_base, params_demand_default_symmetric, params_demand_default_asymmetric_single_cohort, params_demand_default_asymmetric
export params_transition_default_two_cohorts, params_transition_default_three_cohorts, params_transition_default_four_cohorts, stationary_params_awareness_default, stationary_params_default, stationary_params_endogenous_default, stationary_settings_default, set_size, industry_params_default, industry_settings_default, industry_equilibrium, industry_equilibrium_one_step, industry_equilibrium_μ_fixed
export transition_operator_base!, get_transition_operator, ℚf_default!
export π_star, Ψ_hat, v_0, steadystate
export p_i_root, g_i, D_g_i, p_i, π_star_i

export demand_symmetric, demand_asymmetric_single_cohort, demand_asymmetric
export get_profit, get_profit_derivative
export compute_price_equilibrium

end # module
