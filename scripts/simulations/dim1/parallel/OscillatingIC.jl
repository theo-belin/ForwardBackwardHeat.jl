
module OscillatingIC

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

	X_min = -1
	X_max = 1
	T_max = 0.5

	# Space discretisation
	h = 0.1
	X = X_min:h:X_max
	grid = simplexgrid(X)

	# Time discretisation
	k = h/40
	T = k:k:T_max

	x0 = (X_max+X_min)/2

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, _ -> 0)
	
	delta = (P1-P0)/10

	u0_r1, v0_r1, lambda0_r1 = IC_random_1D(
		gamma_0, kappa; 
		x_min = X_min, x_max = X_max, 
		v_min = P0-1, v_max = P1 + 2.5,
		seed = 0
	)
	
	u0_r2, v0_r2, lambda0_r2 = IC_random_1D(
		gamma_0, kappa; 
		x_min = X_min, x_max = X_max, 
		v_min = P0-1, v_max = P1 + 4, 
		seed = 1
	)

	u0_sin, v0_sin, lambda0_sin = IC_sinusoidal_1D(
		gamma_0, kappa; 
		amp = 2*(P1-P0), freq = 4*2*pi, off_set = 1
	)

	tu_r1, tv_r1, tlambda_r1, _ = fbheat(
		u0_r1, v0_r1, lambda0_r1, 
		X, T, problem_data
	)

	tu_r2, tv_r2, tlambda_r2, _ = fbheat(
		u0_r2, v0_r2, lambda0_r2, 
		X, T, problem_data
	)

	tu_sin, tv_sin, tlambda_sin, _ = fbheat(
		u0_sin, v0_sin, lambda0_sin, 
		X, T, problem_data
	)

	if visualize
		p = GridVisualizer(; Plotter = PyPlot, layout = (3, 2), fast = true)
		for i in 1:10:length(tu_r1.t)
			time = tu_r1.t[i]
			scalarplot!(
				p[1, 1], grid, tu_r1[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 1], grid, tv_r1[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[3, 1], grid, tlambda_r1[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[1, 2], grid, tu_sin[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 2], grid, tv_sin[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[3, 2], grid, tlambda_sin[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			reveal(p)
			sleep(1.0e-2)
		end
	end
	if test
		print("OscillatingIC.jl successful")
	end

	tu_r1, tv_r1, tlambda_r1, tu_sin, tv_sin, tlambda_sin
end

end