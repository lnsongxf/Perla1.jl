function solve_transition_dynamics(Q, f_0, T; dt = 0.1)
    # Solve transition dynamics given AwarenessModel object Q

    # generate the corresponding operator for the model
    operator = generate_operator(Q)
    # definte the corresponding ODE problem
    prob = ODEProblem(operator,f_0,(0.0,T))
    # solve the model
    solve(prob, HochOst4(krylov=true), dt=dt)
end