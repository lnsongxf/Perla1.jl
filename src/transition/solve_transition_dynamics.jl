function solve_transition_dynamics(Q!, params, f_0, T; dt = 0.1)
    # Solve transition dynamics given dynamics Q!

    # generate the corresponding operator for the model
    operator = MatrixFreeOperator(Q!, (1, 1)) # extra argument args=(1,1) is dummy
    # definte the corresponding ODE problem
    prob = ODEProblem(operator,f_0,(0.0,T), params)
    # solve the model
    solve(prob, HochOst4(krylov=true), dt=dt)
end