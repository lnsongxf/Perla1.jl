function get_Q(N, μ, θ, θ_d)
    f0(a) = (θ_d + θ) / (θ_d + θ*exp((θ_d + θ)*a)) # (A.4) in appendix
    Q_base = get_Q_base(N, μ, θ)

    # Q(a); note that the only time dependent part is on the first row 
    function Q(a) # (A.3) in appendix
        Q_base[1,1] = -(θ + θ_d*(1-f0(a)))
        Q_base[1,2] = θ + θ_d*(1-f0(a))
        return Q_base
    end
    
    return Q

end

function get_Q_base(N, μ, θ)
    # construct (N+1) by (N+1) base matrix such that
    # the last N rows of base_matrix correspond to the last N rows of Q(a) 
    # based on (A.3) in appendix
    dl = fill(μ, N)
    d = collect(-μ.-(N:-1:0).*(θ/N))
    du = collect((N:-1:1).*(θ/N))
    return SparseArrays.spdiagm(-1 => dl, 0 => d, 1 => du)
end