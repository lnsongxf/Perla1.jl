module Perla1

using SparseArrays, DifferentialEquations, LinearAlgebra, ExponentialUtilities, DiffEqOperators
import Parameters: @with_kw, @unpack

include("transition/model.jl")
include("transition/solve_transition_dynamics.jl")
export AwarenessModel, solve_transition_dynamics

end # module
