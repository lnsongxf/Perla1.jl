function get_Q(N, μ, θ, θ_d)
    f_0(a) = (θ_d + θ) / (θ_d + θ*exp((θ_d + θ)*a)) # (A.4) in appendix
    Q_base = get_Q_base(N, μ, θ)

    # Q(a); note that the only time dependent part is on the first row 
    function Q(a) # (A.3) in appendix
        # FIXIT: there must be a more efficient way to do it
        # deepcopy is needed anyway, since manipulation has to be done
        # if we stick with using Q_base,
        # which is less expensive than computing it every time Q(a) is computed
        Q_a = deepcopy(Q_base)
        Q_a[1,1] = -(θ + θ_d*(1-f_0(a)))
        Q_a[1,2] = θ + θ_d*(1-f_0(a))
        return Q_a
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
    return LinearAlgebra.Tridiagonal(dl,d,du)
end