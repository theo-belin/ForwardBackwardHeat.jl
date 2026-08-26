function heat_equation(u_0, sgrid, T; bc = "Neumann")

    function flux(y, u, edge)
        y[1] = u[1, 1] - u[1, 2]
    end

    function storage(y, u, node)
        y .= u
    end

    function bcondition!(y, u, node)
        if bc == "Dirichlet"
            boundary_dirichlet!(y, u, node; value = 0)
        else
            boundary_neumann!(y, u, node; value = 0)
        end
    end

    sys = VoronoiFVM.System(
        sgrid;
        flux, 
        storage, 
        species = [1], 
        bcondition = bcondition!
        )

    sol = unknowns(sys)
    sol[1,:] .= map(u_0, sgrid)
    tsol = TransientSolution(T[1], sol)
    nT = length(T)
    for iT = 2:nT
        Δt = T[iT] - T[iT-1]
        # solve one timestep via Newton method
        sol = VoronoiFVM.solve(
            sol, sys; 
            time = T[iT], tstep = Δt, 
            called_from_API = true
            )

        append!(tsol, T[iT], sol)
    end
    tsol, sys
end

function nonlinear_diffusion(u_0, phi, sgrid, T)
    function flux(y, u, edge)
        y[1] = phi(u[1, 1]) - phi(u[1, 2])
    end
    function storage(y, u, node)
        y .= u
    end
    sys = VoronoiFVM.System(
        sgrid; 
        flux, 
        storage, 
        species = [1]
        )

    sol = unknowns(sys)
    sol[1,:] .= map(u_0, sgrid)
    tsol = TransientSolution(T[1], sol)
    nT = length(T)
    for iT = 2:nT
        Δt = T[iT] - T[iT-1]
        # solve one time step via Newton method
        sol = VoronoiFVM.solve(
            sol, sys; 
            time = T[iT], tstep = Δt, 
            called_from_API = true
            )
        append!(tsol, T[iT], sol)
    end
    tsol, sys
end

## Solving for a general non-linear phi, in terms of gamma_eps and kappa. 

function fbheat(u0, v0, lambda0, sgrid, T, data)
    """
        u_0, lambda_0: compatible initial conditions.
        phi_eps: (u, lambda) |-> v, the generalized potential (regularized or non-regularized).
        kappa: the non-linear storage term.
        sgrid: the spatial, simplex grid.
        T: the time discretization.

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
        y[1] = data.source(node[1])
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