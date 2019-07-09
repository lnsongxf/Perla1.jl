function E_by_n̂(g, f_a)
   N = length(f_a) - 1
   return dot(f_a[2:end], g.(1:N)) / (1-f_a[1])
end

# solve transition dynamics with transition matrix ℚ
function solve_awareness_dynamics(params, a; tstops = [0; 100.0])
    @unpack N, σ, κ, μ  = params
   Ψ(μ, n) = params.Ψ(μ, n, params)
   f0 = [1.0; zeros(N)]
   tstops = sort(vcat(a, tstops)) # a is included in tstops
   prob = ODEProblem(params.ℚf!,f0,(0.0,tstops[end]),params)
   sol = solve(prob, tstops = tstops)
   f = sol.(a)
   f_0 = (f -> f[1]).(f) # extract f_0 for each a

   # for expectation
   E_n̂ = similar(a)
   E_n̂_m1 = similar(a)
   E_n̂_m2 = similar(a)
   for a_index in eachindex(a)
      f_a = f[a_index]
      E_n̂[a_index] = E_by_n̂(n -> Ψ(μ, n)^(σ*(κ-1)), f_a)
      E_n̂_m1[a_index] = E_by_n̂(n -> n*Ψ(μ, n)^(σ*(κ-1)-1), f_a)
      E_n̂_m2[a_index] = E_by_n̂(n -> n*Ψ(μ, n)^(σ*(κ-1)-2), f_a)
   end
   return (E_n̂ = E_n̂, E_n̂_m1 = E_n̂_m1, E_n̂_m2 = E_n̂_m2, f_0 = f_0, f = f)
end

# (A.7)
q(E_n̂) = E_n̂ 

# (A.8)
function Upsilon(E_n̂_m1, E_n̂_m2, params) 
   @unpack μ, σ, κ = params
   return 1 + σ / (1 - (1 - σ * (κ - 1)) * E_n̂_m2/E_n̂_m1)
end

# (A.9)
function Q(E_n̂, E_n̂_m1, E_n̂_m2, f_0, E, params) 
   @unpack μ, σ, κ  = params
   return (E*((1 .- f_0) .* Upsilon.(E_n̂_m1, E_n̂_m2, Ref(params)).^(1-κ).*q.(E_n̂)))^(1/(κ-1))
end

# (A.10)
function B(E_n̂, E_n̂_m1, E_n̂_m2, f_0, E, params) 
   @unpack μ, σ, κ  = params
   return (E*((1 .- f_0) .* Upsilon.(E_n̂_m1, E_n̂_m2, Ref(params)).^(-κ).*q.(E_n̂))) /
   (E*((1 .- f_0) .* Upsilon.(E_n̂_m1, E_n̂_m2, Ref(params)).^(1-κ).*q.(E_n̂)))
end

# (B.7)
# value of the firm choosing θ_i when all other firms choose θ
 function π_star(a, p, f_0, E_n̂, mc, Y, params)
   @unpack σ, κ, N = params
   return Y*(p - mc)*(1 - f_0)/N*p^(-κ) * q(E_n̂) 
end

# (A.22)
function v_0(a, p, f_0, E_n̂, mc, Y, pdv_weights, params, settings)
   @unpack δ_M, ρ, N = params
   return pdv_weights ⋅ π_star.(a, p, f_0, E_n̂, mc, Y, Ref(params))
end

function steadystate(params, settings)
   @unpack σ, κ, ρ, α, N, δ_M, δ_k, θ, μ, z_M = params

   # age distribution and expectation setups
   age_dist = Exponential(1/δ_M)
   E_age_dist = expectation(age_dist, n = settings.expectation_node_count)

   # solve awareness dynamics
   a = nodes(E_age_dist) # extract a

   @unpack E_n̂, E_n̂_m1, E_n̂_m2, f_0, f = solve_awareness_dynamics(params, a; tstops = a);

   Q_val = Q(E_n̂, E_n̂_m1, E_n̂_m2, f_0, E_age_dist, params);
   B_val = B(E_n̂, E_n̂_m1, E_n̂_m2, f_0, E_age_dist, params);

   # solve k, M
   function resid!(resid, k_and_M)
      k = k_and_M[1]
      M = k_and_M[2]
      if (k < 0 || M < 0)
         resid[:] = [Inf; Inf ]
         return
      end

      # (A.6) and (A.7)
      resid[:] = [δ_M - δ_k - Q_val*B_val^(-1)*(k^α)*M^(1/(κ-1))*(z_M/(M*(κ-1)) - α/k);
                  ρ + δ_k - α*Q_val*B_val^(-1)*M^(1/(κ-1))*k^(α-1)]
   end
   solved = nlsolve(resid!, settings.ss_iv(params), autodiff = :forward)
   k, M = solved.zero

   # post calculation Q_valuantities, (A.13) - (A.21)
   C = Q_val/B_val*M^(1/(κ-1))*k^α - δ_k*k - δ_M*M/z_M # (A.13)
   capital_share = α*B_val # (A.14)
   labor_share = (1-α)*B_val # (A.15)
   profit_share = 1 - B_val # (A.16)
   mc = M^(1/(κ-1))*Q_val # (A.17)
   Z = M^(1/(κ-1))*Q_val/B_val # (A.18)
   w = (1-α)*Z*B_val*k^α # (A.19)
   Y = Z*k^α # (A.20)
   r = ρ + δ_M
   TobinsQ = 1 + (1-B_val)/(1 - r)*Y/k # (A.21)
   p = Upsilon.(E_n̂_m1, E_n̂_m2, Ref(params)) .* mc # (A.8) symmetric prices at each a
   pdv_weights = weights(E_age_dist) .* exp.((δ_M-r) * a) / δ_M # for calculating int e^(-r a) f(a) da using existing quadrature nodes

   return (k = k, M = M, Q = Q_val, B = B_val, C = C, capital_share = capital_share, labor_share = labor_share, profit_share = profit_share, mc = mc, Z = Z, w = w, Y = Y, TobinsQ = TobinsQ,
   E_n̂ = E_n̂, E_n̂_m1 = E_n̂_m1, E_n̂_m2 = E_n̂_m2,
           f_0 = f_0, f = f, p = p, a = a,
           q = q.(E_n̂),
           v_0 = v_0(a, p, f_0, E_n̂, mc, Y, pdv_weights, params, settings), E_age_dist = E_age_dist, pdv_weights = pdv_weights,
           )
end
