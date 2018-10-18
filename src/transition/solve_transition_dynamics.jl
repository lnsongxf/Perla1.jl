function solve_transition_dynamics(Q, f_0)
    # Solve transition dynamics given AwarenessModel object Q
    # using Krylov methods
    function sol_krylov(t)
        Q.t = t
        return expv(t,Q,f_0)
    end
    return sol_krylov;
end