## Solver for piece-wise linear phi

function fbheat(u_0, psi_0, phi0, phi1, P0, P1, sgrid, T)

    function phi_f(u, psi)
        (1 - psi) * phi0(u) + psi * phi1(u)
    end

    function psi_f(u, psi)
        v = phi_f(u, psi)
        if v > P1
            1.0
        elseif v < P0
            0.0
        else
            psi
        end
    end

    function flux(y, u, edge)
        K = edge.node[1]
        L = edge.node[2]
        y[1] = phi_f(u[1, 1], psi[K]) - phi_f(u[1, 2], psi[L])
    end
    
    function storage(y, u, node)
        y .= u
    end

    sys = VoronoiFVM.System(sgrid; flux, storage, species = [1])
    sol = unknowns(sys)
    psi = unknowns(sys)
    sol[1, :] .= map(u_0, sgrid)
    psi[1, :] = map(psi_0, sgrid)
    tu = TransientSolution(T[1], sol)
    tpsi = TransientSolution(T[1], psi)
    nT = length(T)
    for iT = 2:nT
        Δt = T[iT] - T[iT-1]
        # solve one timestep via Newton method
        control = VoronoiFVM.SolverControl()
        #control.verbose=""
        sol = VoronoiFVM.solve(
            sol,
            sys;
            time = T[iT],
            tstep = Δt,
            called_from_API = true,
            control,
        )
        psi[1,:] = psi_f.(sol[1,:], psi[1,:])
        append!(tu, T[iT], sol)
        append!(tpsi, T[iT], copy(psi))
    end
    tu, tpsi, sys
end


function fbheat_plotnikov(u_0, lambda_0, eps, P0, P1, alpha, kappa, kappa_0, sgrid, T)

    function gamma_eps(v, lambda)
        if v <= P0 - eps
            0
        elseif P0 - eps < v && v < P0
            lambda/eps*(v-P0) + lambda
        elseif P0 <= v && v <= P1
            lambda
        elseif P1 < v && v < P1 + eps
            (1-lambda)/eps*(v-P1) + lambda
        else
            1
        end
    end

    #Inverse function of v |--> alpha*v + beta*gamma_eps(v,lambda) + kappa_0
    function phi_eps(u, lambda)
        if u <= alpha*(P0 - eps) + beta_0
            (u - beta_0)/alpha
        elseif alpha*(P0-eps) + beta_0 < u && u < alpha*P0 + beta*lambda + beta_0
            (beta*lambda*P0 + eps*(u -beta*lambda - beta_0))/(beta*lambda + eps*alpha)
        elseif alpha*P0 + beta*lambda + beta_0 <= u && u <= alpha*P1 + beta*lambda + beta_0
            (u - beta*lambda - beta_0)/alpha
        elseif alpha*P1 + beta*lambda + beta_0 < u && u < alpha*P1 + beta + beta_0
            (beta*(1-lambda)*P1 + eps*(u - beta*lambda - beta_0))/(beta*(1-lambda) + eps*alpha)
        else
            (u - beta - beta_0)/alpha
        end
    end

    function flux(y, u, edge)
        K = edge.node[1]
        L = edge.node[2]
        y[1] = phi_eps(u[1, 1], lambda[K]) - phi_eps(u[1,2], lambda[L])
    end
    
    function storage(y, u, node)
        y .= u
    end

    sys = VoronoiFVM.System(sgrid; flux, storage, species = [1])
    sol = unknowns(sys)
    sol[1, :] .= map(u_0, sgrid)
    lambda::Vector{Float64} = map(lambda_0, sgrid)
    tu = TransientSolution(T[1], sol)
    tlambda = TransientSolution(T[1], lambda)
    nT = length(T)
    for iT = 2:nT
        Δt = T[iT] - T[iT-1]
        # solve one timestep via Newton method
        control = VoronoiFVM.SolverControl()
        #control.verbose=""
        sol = VoronoiFVM.solve(
            sol,
            sys;
            time = T[iT],
            tstep = Δt,
            called_from_API = true,
            control,
        )
        v = phi_eps.(sol[1,:], lambda)
        lambda = (sol[1, :] .- (alpha.*v .+ beta_0))./beta
        append!(tu, T[iT], sol)
        append!(tlambda, T[iT], copy(lambda)) 
    end
    tu, tlambda, sys
end


