# Solving for a general non-linear phi, in terms of gamma_eps and kappa. 

function fbheat(u0, v0, lambda0, sgrid, T, data)
    """
        u_0, v0, lambda_0: compatible initial conditions.
        sgrid: the spatial, simplex grid.
        T: the time discretization.
        data: struct containing 
            gamma_eps: (v,lambda) |-> lambda, the map for the hysteresis operator, 
            kappa: (v, lambda) |-> u, the non-linear storage term,
            phi_eps: (u, lambda) |-> v, the relaxed potential (regularized or non-regularized).

            Finite volume solver for the entropy solution of the forward backward nonlinear heat equation. 
    """
    gamma_eps, kappa, phi_eps = data.gamma_eps, data.kappa, data.phi_eps
    kappa0(v) = kappa(v,0.0)
    kappa1(v) = kappa(v,1.0)

    function flux!(y, u, edge, data)
        K = edge.node[1]
        L = edge.node[2]
        y[1] = phi_eps(u[1, 1], lambda[K]) - phi_eps(u[1,2], lambda[L])
        return nothing
    end
    
    function storage!(y, u, node, data)
        y .= u
        return nothing
    end

    function source!(y, node, data)
        y[1] = data.source(node)
    end

    physics = VoronoiFVM.Physics(;
        flux = flux!,
        storage = storage!, 
        source = source!,
        data = data
    )

    sys = VoronoiFVM.System(sgrid, physics; assembly = :edgewise)

    enable_species!(sys, 1, [1])

    u = unknowns(sys)
    v = unknowns(sys)
    lambda = unknowns(sys)
    
    u[1,:] = map(u0, sgrid)
    v[1,:] = map(v0, sgrid)
    lambda[1,:] = map(lambda0, sgrid)
    
    tu = TransientSolution(T[1], u)
    tv = TransientSolution(T[1], v)
    tlambda = TransientSolution(T[1], lambda)
   
    nT = length(T)
    for iT = 2:nT
        Δt = T[iT] - T[iT-1]
        # solve one time step via Newton method
        control = VoronoiFVM.SolverControl(maxiters = 200)
        u = VoronoiFVM.solve(
            sys;
            inival = u,
            time = T[iT],
            tstep = Δt,
            called_from_API = true,
            control
        )
        append!(tu, T[iT], u)
        v[1,:] = phi_eps.(u[1,:], lambda[1,:])
        append!(tv, T[iT], copy(v))
        #lambda[1,:] = gamma_eps.(v[1,:], lambda[1,:])
        lambda[1,:] = (u[1,:] .- kappa0.(v[1,:]))./(kappa1.(v[1,:]) .- kappa0.(v[1,:]))
        append!(tlambda, T[iT], copy(lambda))
    end
    tu, tv, tlambda, sys
end