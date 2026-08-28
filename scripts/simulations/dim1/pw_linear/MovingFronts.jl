module MovingFronts

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

function main(; visualize = true, test = false)

	X_min = -1
	X_max = 1
	T_max = 0.5

	# Space discretization
	h = 0.1
	X = X_min:h:X_max
	sgrid = simplexgrid(X)

	# Time discretization
	k = h/80
	T = k:k:T_max

	x0 = (X_max+X_min)/2

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, _ -> 0)
	
	# Riemann data with saturated phase fraction
	u0_sat, v0_sat, lambda0_sat = IC_riemann_1D(
		gamma_0, kappa,
		P1 + 1.9, P1, 
		1, 0; 
		x0 = x0
	)

	# Solver
	tu_sat, tv_sat, tlambda_sat, _ = fbheat(
		u0_sat, v0_sat, lambda0_sat, 
		sgrid, T, problem_data
	)

	# Riemann data with non-saturated phase fraction
	u0_unsat, v0_unsat, lambda0_unsat = IC_riemann_1D(
		gamma_0, kappa, 
		P1 + 1.9, P1, 
		1, 0.5; 
		x0 = x0
	)
	
	# Solver
	tu_unsat, tv_unsat, tlambda_unsat, _ = fbheat(
		u0_unsat, v0_unsat, lambda0_unsat, 
		sgrid, T, problem_data
	)

	print(u0_sat.(X))
	print(v0_sat.(X))
	print(lambda0_sat.(X))

	print(tu_sat[1,:,1])
	print(tv_sat[1,:,1])
	print(tlambda_sat[1,:,1])

	# Step-by-step visualization
	if visualize
		p = GridVisualizer(; Plotter = PyPlot, layout = (2, 3), fast = true)
		n_t = length(tu_sat.t)
		log_times = Int.(floor.((1.5).^(1:1:(log(n_t)/log(1.5)))))
		for i = log_times
			time = tu_sat.t[i]
			scalarplot!(
				p[1, 1], sgrid, tu_sat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[1, 2], sgrid, tv_sat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[1, 3], sgrid, tlambda_sat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :blue, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 1], sgrid, tu_unsat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :red, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 2], sgrid, tv_unsat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :red, label = "numerical",
				markershape = :circle, markevery = 1
			)
			scalarplot!(
				p[2, 3], sgrid, tlambda_unsat[1, :, i]; title = @sprintf("t=%.3g", time),
				color = :red, label = "numerical",
				markershape = :circle, markevery = 1
			)
			reveal(p)
			sleep(1.0)
		end
	end

	if test
		print("MovingFronts.jl successful")
	end

	tu_sat, tv_sat, tlambda_sat, tu_unsat, tv_unsat, tlambda_unsat
end

end