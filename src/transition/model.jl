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
    # transition generator info
    # ==================================
    dl # lower diagonal elements of Q
    d # diagonal elements of Q
    du # upper diagonal elements of Q

    # ==================================
    # time dependency
    # ==================================
    t::Float64 # time
    # ==================================
    # constructor
    # ==================================
    function AwarenessModel(N, μ, θ, θ_d, f0)
        dl, d, du = get_Q_base(N, μ, θ)
        new(N, μ, θ, θ_d, f0, dl, d, du, 0.0)
    end
    function AwarenessModel(N, μ, θ, θ_d, f0, t)
        dl, d, du = get_Q_base(N, μ, θ)
        new(N, μ, θ, θ_d, f0, dl, d, du, t)
    end
end

# overload size, ishermitian, mul!, opnorm from LinearAlgebra to use expv
LinearAlgebra.size(Q::AwarenessModel, i::Int64) = (Q.N+1)
LinearAlgebra.ishermitian(Q::AwarenessModel) = false

function LinearAlgebra.mul!(y, Q::AwarenessModel, b)
    # note that the operator is transposed
    N = Q.N
    y[1] = -(Q.θ + Q.θ_d*(1-Q.f0(Q.t)))*b[1]+Q.dl[1]*b[2]
    for i in 2:N
        y[i] = Q.du[i-1]*b[i-1] + Q.d[i]*b[i] + Q.dl[i]*b[i+1]
    end
    y[2] = (Q.θ + Q.θ_d*(1-Q.f0(Q.t)))*b[1]+Q.d[2]*b[2]+Q.dl[2]*b[3]
    y[end] = Q.du[end]*b[N] + Q.d[end]*b[N+1]
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