# solve transition dynamics with transition operator O! (either matrix or matrix-free)
function solve_transition_dynamics(O!, params, f0, T; dt = 0.25)
    @unpack θ, cohorts, ts_cohort, ζ = params
    @assert length(unique(cohorts)) == 1 # assert that cohort sizes are invariant over t.
    N = cohorts[1]
    prob = ODEProblem(O!,f0,(0.0,T), params)
    solve(prob, LinearExponential(krylov=:simple), tstops=union(0:dt:T, ts_cohort, T))
end

# solve transition dynamics with transition matrix Q 
function solve_transition_dynamics_single_cohort(Q, f0, T; dt = 0.25)
    df(f,p,a) = Q(a, f)' * f
    prob = ODEProblem(df,f0,(0.0,T))
    solve(prob)
end
