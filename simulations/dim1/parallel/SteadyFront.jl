module SteadyFront

using Printf

using ForwardBackwardHeat
using ForwardBackwardHeat.Parallel
using ExtendableGrids, Plots
using GridVisualize
using VoronoiFVM

## Problem data structure to avoid global variables
mutable struct ProblemData
    phi::Any
	kappa::Any
	gamma_eps::Any
	phi_eps::Any
	source::Any
end

function main(; visualize = true, test = false)

	X_min = -1
	X_max = 1
	T_max = 0.5

	# Space discretisation
	h = 0.2
	X = X_min:h:X_max
	sgrid = simplexgrid(X)

	# Time discretisation
	k = h/80
	T = k:k:T_max

	x0 = (X_max+X_min)/2

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, _ -> 0)
	
	delta = (P1-P0)/10
	u0, v0, lambda0 = IC_riemann_1D(
		gamma_0, kappa,
		P1 - delta, P0 + delta, 
		0, 1; 
		x0 = x0
	)

	tu, tv, tlambda, _ = fbheat(
		u0, v0, lambda0, 
		sgrid, T, problem_data
	)

	if visualize
		p = GridVisualizer(; Plotter = Plots, layout = (3, 1), fast = true)
		for i in 1:10:length(tu.t)
			time = tu.t[i]
			scalarplot!(
				p[1, 1], sgrid, tu[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerica",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 1], sgrid, tv[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerica",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[3, 1], sgrid, tlambda[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerica",
				markershape = :circle, markevery = 1
			)
			reveal(p)
			sleep(1.0)
		end
    end
	if test
		print("SteadyFront.jl successfu")
	end
	
	tu, tv, tlambda

end

end