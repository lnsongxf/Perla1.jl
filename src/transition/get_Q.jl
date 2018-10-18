function get_Q_base(N, μ, θ)
    # construct (N+1) by (N+1) base matrix such that
    # the last N rows of base_matrix correspond to the last N rows of Q(a) 
    # based on (A.3) in appendix
    dl = fill(μ, N)
    d = collect(-μ.-(N:-1:0).*(θ/N))
    du = collect((N:-1:1).*(θ/N))
    return dl, d, du
end