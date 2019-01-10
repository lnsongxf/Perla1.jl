module Perla1

using SparseArrays, DifferentialEquations, LinearAlgebra, ExponentialUtilities, DiffEqOperators, Parameters, SpecialFunctions
import Parameters: @with_kw, @unpack

include("transition/solve_transition_dynamics.jl")
include("demand/params.jl")
export solve_transition_dynamics
export params_demand_default_base, params_demand_default_symmetric, params_demand_default_asymmetric_single_cohort, params_demand_default_asymmetric

end # module
