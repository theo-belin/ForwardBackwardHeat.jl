module ConvergenceError

using Printf
using PyPlot
using ForwardBackwardHeatProject
using ForwardBackwardHeatProject.PWLinear
using ExtendableGrids
using GridVisualize
using VoronoiFVM

mutable struct ProblemData
    phi::Any
	kappa::Any
	gamma_eps::Any
	phi_eps::Any
	source::Any
end

function main(; visualize = true, test = false, p = 2)

	# Space-time domain
	X_min = -1
	X_max = 1

	x0 = (X_max+X_min)/2

	T_max = 0.5

	function source(x)
		pi^2*cos(pi*x)
	end

	# Create problem data structure
	problem_data = ProblemData(phi, kappa, gamma_0, phi_0, source)

	# Steady sate with the above source term 
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

	u0, v0, lambda0 = prepare_IC(
		gamma_0, kappa, 
		v0_ini, lambda0_ini
	)
	
	# Sequence of mesh sizes
	H = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
	error = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	# Exponent of Lebesgue space
	p = 2

	for i = 1:length(H)
		h = H[i]

		# Space discretisation
		X = X_min:h:X_max
		grid = simplexgrid(X)

		# Time discretisation
		k = h/40
		T = k:k:T_max

		problem_data = ProblemData(phi, kappa, gamma_0, phi_0, source)

		tu, tv, tlambda, sys = fbheat(
			u0, v0, lambda0, 
			X, T, problem_data
		)

		er = 0
		for j in 1:(length(tu.t)-1)
			tj = tu.t[j]
			tjp = tu.t[j + 1]
			dt = tjp - tj
			er = er + dt*VoronoiFVM.lpnorm(sys, tu(tj) - tu(k), p)^p
		end
		error[i] = (er)^(1/p)
	end

	if visualize
		p = GridVisualizer(; Plotter = PyPlot, layout = (1, 1), fast = true)
		scalarplot!(
			p, log10.(H.^(-1)), log10.(error.^(-1));
			title = L"Convergence error in $L^p$", size = (250, 100), 
			markershape = :circle, markevery = 1
			)
		reveal(p)
	end

	if test
		print("ConvergenceError.jl successful")
	end

	H, error

end

end