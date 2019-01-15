using Test
using Perla1
using LinearAlgebra
using DifferentialEquations
using DiffEqOperators
import Parameters: @with_kw, @unpack

include("demand/runtests.jl")
include("transition/runtests.jl")
