# solve transition dynamics with transition operator O! (either matrix or matrix-free)
function solve_transition_dynamics(O!, params, f_0, T; dt = 0.25)
    # unpack params
    @unpack θ, cohorts, ts_cohort, μ = params
    N = cohorts[1] # assume that N_t is invariant across all t
        
    # definte the corresponding ODE problem
    prob = ODEProblem(O!,f_0,(0.0,T), params)
    
    # solve the model
    solve(prob, LinearExponential(krylov=:simple), tstops=union(0:dt:T, ts_cohort, T))
end