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
    # time dependency
    # ==================================
    t::Float64 # time
    # ==================================
    # constructor
    # ==================================
    function AwarenessModel(N, μ, θ, θ_d, f0)
        new(N, μ, θ, θ_d, f0, 0.0)
    end
    function AwarenessModel(N, μ, θ, θ_d, f0, t)
        new(N, μ, θ, θ_d, f0, t)
    end
end

# overload size, ishermitian, mul!, opnorm from LinearAlgebra to use expv
LinearAlgebra.size(Q::AwarenessModel, i::Int64) = (Q.N+1)
LinearAlgebra.ishermitian(Q::AwarenessModel) = false

function LinearAlgebra.mul!(y, Q::AwarenessModel, b)
    # params = (N = N, μ = μ, θ = θ, θ_d = θ_d, f0 = f0) # TODO: later add this to unit tests
    params = (N = Q.N, μ = Q.μ, θ = Q.θ, θ_d = Q.θ_d, f0 = Q.f0) # TODO: later replace this with params
    @unpack N, μ, θ, θ_d, f0 = params
    y[1] = -(θ + θ_d*(1-f0(Q.t)))*b[1] + μ*b[2]
    for i in 2:N
        y[i] = θ*((N+2-i)/N)*b[i-1] - (μ+θ*((N+1-i)/N))*b[i] + μ*b[i+1]
    end
    y[2] = (θ + θ_d*(1-f0(Q.t)))*b[1] - (μ+θ*((N-1)/N))*b[2] + μ*b[3]
    y[end] = (θ/N)*b[N] - μ*b[N+1]
end

# for multiple cols
function LinearAlgebra.mul!(Y, Q::AwarenessModel, B::AbstractMatrix{T}) where {T}
    N = Q.N
    for j in 1:(N+1)
        y = Y[:,j]
        mul!(y, Q, B[:,j])
        Y[:,j] = y
    end
end

function LinearAlgebra.opnorm(Q::AwarenessModel, p)
    # NOTE: ASSUMES μ,θ > 0 AND p = Inf
    # Inf norm of Q' is equivalent to the max_i (sum of ith row in abs vals) of Q'
    # Except the very first two rows, this is zero for every row; see (A.3) in appendix 
    # Hence, sufficient to consider max of the first two rows
    # [[-(θ + θ_d*(1-f0)) μ];
    # [(θ + θ_d*(1-f0)) -μ-(N-1)*θ/N]]
    ii = abs(Q.θ + Q.θ_d*(1-Q.f0(Q.t)))
    return (max(ii + Q.μ, ii + Q.μ + (Q.N-1)*Q.θ/Q.N))
end

# return corresponding operator for ODE solvers
function generate_operator(Q::AwarenessModel) 
    function awareness_operator_basis(df,f,p,t)
        Q.t = t
        mul!(df, Q, f)
    end
    operator = MatrixFreeOperator(awareness_operator_basis, (1,1)) # extra argument p=(1,1) is dummy
    return operator
end
DiffEqBase.isinplace(::MatrixFreeOperator, num) = true