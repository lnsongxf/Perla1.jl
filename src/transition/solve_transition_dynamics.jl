function solve_transition_dynamics(Q, f_0, T)
    # solve transition dynamics given 
    # Q; N by N matrix generator
    # f_0; N vector of initial distribution
    # T; Float64 terminal time
    df(f,p,a) = Q(a)' * f
    prob = ODEProblem(df,f_0,(0.0,T))
    return solve(prob);
end