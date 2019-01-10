module Perla1

using SparseArrays, DifferentialEquations, LinearAlgebra, ExponentialUtilities, DiffEqOperators, Parameters, SpecialFunctions, OffsetArrays
using ForwardDiff: jacobian, derivative
import Parameters: @with_kw, @unpack

include("transition/solve_transition_dynamics.jl")
include("demand/params.jl")
include("demand/demand.jl")
include("demand/profit.jl")

export solve_transition_dynamics
export params_demand_default_base, params_demand_default_symmetric, params_demand_default_asymmetric_single_cohort, params_demand_default_asymmetric
export demand_symmetric, demand_asymmetric_single_cohort, demand_asymmetric
export get_profit, get_profit_derivative

end # module
