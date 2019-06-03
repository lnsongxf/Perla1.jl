# solve transition dynamics with transition operator O! (either matrix or matrix-free)
function solve_transition_dynamics(O!, params, f_0, T; dt = 0.25)
    @unpack θ, cohorts, ts_cohort, μ = params
    @assert length(unique(cohorts)) == 1 # assert that cohort sizes are invariant over t.
    N = cohorts[1]
    prob = ODEProblem(O!,f_0,(0.0,T), params)
    solve(prob, LinearExponential(krylov=:simple), tstops=union(0:dt:T, ts_cohort, T))
end
