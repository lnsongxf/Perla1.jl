module Perla1

using SparseArrays, DifferentialEquations, LinearAlgebra, ExponentialUtilities

include("transition/model.jl")
include("transition/get_Q.jl")
include("transition/solve_transition_dynamics.jl")
export AwarenessModel, solve_transition_dynamics

end # module
