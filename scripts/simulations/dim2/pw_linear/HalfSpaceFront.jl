module HalfSpaceFront

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

	# Space discretization
	h = 0.1
	X = -1:h:1
    Y = -1:h:1

    sgrid = simplexgrid(X, Y)

	# Time discretization
	k = h/10
	T= 0.1:k:1

	chi_halfspace = construct_2D_indicator_fun("Half-plane")

	## Create problem data structure
    problem_data = ProblemData(phi, kappa, gamma_0, phi_0, _ -> 0)
	
	# Riemann data with saturated phase fraction in the disk
	u0_sat, v0_sat, lambda0_sat = IC_riemann_2D(
		gamma_0, kappa,
		P1 + 1.9, P1, 
		1, 0,
		chi_halfspace
	)

	# Solver
	tu_sat, tv_sat, tlambda_sat, _ = fbheat(
		u0_sat, v0_sat, lambda0_sat, 
		sgrid, T, problem_data
	)

	# Riemann data with non-saturated phase fraction
	u0_unsat, v0_unsat, lambda0_unsat = IC_riemann_2D(
		gamma_0, kappa, 
		P1 + 1.9, P1, 
		1, 0.5,
		chi_halfspace
	)
	
	# Solver
	tu_unsat, tv_unsat, tlambda_unsat, _ = fbheat(
		u0_unsat, v0_unsat, lambda0_unsat, 
		sgrid, T, problem_data
	)

	# Step-by-step visualization
	if visualize
		p = GridVisualizer(; Plotter = PyPlot, layout = (2, 3), legend =:rt, fast = true, size = (1500, 600), title = "Expanding Disk with Saturated (top) and Unsaturated (bottom) Phase")
		n_t = length(tu_sat.t)
		log_times = Int.(floor.((1.5).^(1:1:(log(n_t)/log(1.5)))))
		for i = log_times
			time = tu_sat.t[i]
			scalarplot!(
				p[1, 1], sgrid, tu_sat[1, :, i]; 
				flimits = (-1.1, 5),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = L"$u$",
			)
			scalarplot!(
				p[1, 2], sgrid, tv_sat[1, :, i]; 
				flimits = (P1, P1 + 1.9),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = L"$v$"
			)
			scalarplot!(
				p[1, 3], sgrid, tlambda_sat[1, :, i];
				flimits = (0, 1), 
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = L"$\lambda$"
			)
			scalarplot!(
				p[2, 1], sgrid, tu_unsat[1, :, i];
				flimits = (-1.1, 5),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = L"$u$"
			)
			scalarplot!(
				p[2, 2], sgrid, tv_unsat[1, :, i]; 
				flimits = (P1, P1 + 1.9),
				levels = 0,
				title = @sprintf("t=%.3g", time), 
				label = L"$v$"
			)
			scalarplot!(
				p[2, 3], sgrid, tlambda_unsat[1, :, i];
				flimits = (0, 1),
				levels = 0,
				title = @sprintf("t=%.3g", time),
				label = L"$\lambda$"
			)
			reveal(p)
			sleep(1.0)
		end
	end

	if test
		print("HalfSpaceFront.jl successful")
	end
end

end