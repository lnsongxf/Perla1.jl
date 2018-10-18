function solve_transition_dynamics(Q, f_0)
    # Solve transition dynamics given AwarenessModel object Q
    # using Krylov methods
    sol_krylov(a) = expv(a,Q,f_0)
    return sol_krylov;
end