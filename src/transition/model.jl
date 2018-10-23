mutable struct AwarenessModel 
    # ==================================
    # params for awareness process
    # ==================================
    N # number of products
    μ
    θ
    θ_d    
    f0::Function
    # ==================================
    # constructor
    # ==================================
    function AwarenessModel(N, μ, θ, θ_d, f0)
        new(N, μ, θ, θ_d, f0)
    end
end

function LinearAlgebra.mul!(y, Q::AwarenessModel, b, t)
    # params = (N = N, μ = μ, θ = θ, θ_d = θ_d, f0 = f0) # TODO: later add this to unit tests
    params = (N = Q.N, μ = Q.μ, θ = Q.θ, θ_d = Q.θ_d, f0 = Q.f0) # TODO: later replace this with params
    @unpack N, μ, θ, θ_d, f0 = params
    y[1] = -(θ + θ_d*(1-f0(t)))*b[1] + μ*b[2]
    for i in 2:N
        y[i] = θ*((N+2-i)/N)*b[i-1] - (μ+θ*((N+1-i)/N))*b[i] + μ*b[i+1]
    end
    y[2] = (θ + θ_d*(1-f0(t)))*b[1] - (μ+θ*((N-1)/N))*b[2] + μ*b[3]
    y[end] = (θ/N)*b[N] - μ*b[N+1]
end

# return corresponding operator for ODE solvers
function generate_operator(Q::AwarenessModel) 
    function awareness_operator_basis(df,f,p,t)
        mul!(df, Q, f, t)
    end
    operator = MatrixFreeOperator(awareness_operator_basis, (1,1)) # extra argument p=(1,1) is dummy
    return operator
end
DiffEqBase.isinplace(::MatrixFreeOperator, num) = true