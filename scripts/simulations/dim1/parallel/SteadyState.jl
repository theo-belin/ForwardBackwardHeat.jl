module SteadyState

using Printf
using PyPlot
using ForwardBackwardHeat
using ForwardBackwardHeat.Parallel
using ExtendableGrids
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

	# Space discretization
	h = 0.1
	X = -1:h:1
    sgrid = simplexgrid(X)

	# Time discretization
	k = h/10
	T= 0.1:k:1

	# Source term
	function source(node)
		x = node[1]
		pi^2*cos(pi*x)
	end

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, source)

	# Steady sate with the source term 
	function v0_ini(x)
		cos(pi*x)
	end

	function lambda0_ini(x)
		if x > -1/2 && x < 1/2
			1
		else
			0
		end
	end

	# preparing IC
	u0, v0, lambda0 = prepare_IC(
		gamma_0, kappa, 
		v0_ini, lambda0_ini
	)

	# Solver
	tu, tv, tlambda, _ = fbheat(
		u0, v0, lambda0, 
		sgrid, T, problem_data
	)
	
	# Step-by-step visualization
	if visualize
		p = GridVisualizer(; Plotter = PyPlot, layout = (3, 1), fast = true)
		n_t = length(tu.t)
		log_times = Int.(floor.((1.5).^(1:1:(log(n_t)/log(1.5)))))
		for i = log_times
			time = tu.t[i]
			scalarplot!(
				p[1, 1], sgrid, tu[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 1], sgrid, tv[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[3, 1], sgrid, tlambda[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			reveal(p)
			sleep(1.0)
		end
	end

	if test
		print("SteadyState.jl successful")
	end
end

end