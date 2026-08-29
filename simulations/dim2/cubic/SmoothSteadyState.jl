module SteadyState

using Printf

using ForwardBackwardHeat
using ForwardBackwardHeat.Cubic
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

	# Space discretization
	h = 0.1
	X = -1:h:1
    Y = -1:h:1

    sgrid = simplexgrid(X, Y)

	# Time discretization
	k = h/10
	T= 0.1:k:1

	# Source term
	function source(node)
		x, y = node[1], node[2]
		2*pi^2*cos(pi*x)*cos(pi*y)
	end

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, source)

	# Steady sate with the source term 
	function v0_ini(x, y)
		cos(pi*x)*cos(pi*y)
	end

	function lambda0_ini(x, y)
		if x > -1/2 && x < 1/2 && y > -1/2 && y < 1/2
			1
		else
			0
		end
	end

	# preparing IC
	u0, v0, lambda0 = prepare_IC_2D(
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
		p = GridVisualizer(; Plotter = Plots, layout = (1, 3), legend =:rt, fast = true, size = (1500, 600), 
		title = "Expanding Disk with Saturated (top) and Unsaturated (bottom) Phase")
		n_t = length(tu.t)
		log_times = Int.(floor.((1.5).^(1:1:(log(n_t)/log(1.5)))))
		for i = log_times
			time = tu.t[i]
			scalarplot!(
				p[1, 1], sgrid, tu[1, :, i]; 
				flimits = (-4, 4),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = "u",
			)
			scalarplot!(
				p[1, 2], sgrid, tv[1, :, i] .- tv[1, :, 1]; 
				flimits = (-1, 1),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = "v"
			)
			scalarplot!(
				p[1, 3], sgrid, tlambda[1, :, i];
				flimits = (0, 1), 
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = "lambda"
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