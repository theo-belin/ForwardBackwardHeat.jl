module MovingFrontsForcing

using Printf
using PyPlot
using ForwardBackwardHeat
using ForwardBackwardHeat.PWLinear
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

function main(; visualize = true)

	X_min = -1
	X_max = 1
	T_max = 1

	x0 = (X_max+X_min)/2

	function source(x)
		if x < x0
			4
		else
			0
		end
	end

	# Create problem data structure
	problem_data = ProblemData(phi, kappa, gamma_0, phi_0, source)

	# Space discretisation
	h = 0.2
	X = X_min:h:X_max
	grid = simplexgrid(X)

	# Time discretisation
	k = h/80
	T = k:k:T_max

	u0, v0, lambda0 = IC_riemann_1D(
		gamma_0, kappa, 
		P1, P1, 
		1, 0; 
		x0 = x0
	)

	tu, tv, tlambda, _ = fbheat(
		u0, v0, lambda0, 
		X, T, problem_data
	)

	if visualize

		p = GridVisualizer(; Plotter = PyPlot, layout = (3, 1), fast = true)
		for i in 1:10:length(tu.t)
			time = tu.t[i]
			scalarplot!(
				p[1, 1], grid, tu[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 1], grid, tv[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[3, 1], grid, tlambda[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			reveal(p)
			sleep(1.0e-2)
		end
    end

	tu, tv, tlambda
end

end