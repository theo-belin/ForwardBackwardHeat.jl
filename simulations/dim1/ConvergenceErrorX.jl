module ConvergenceErrorX

using Printf

using ForwardBackwardHeat
using ExtendableGrids, Plots
using GridVisualize
using VoronoiFVM

mutable struct ProblemData
    phi::Any
	kappa::Any
	gamma_eps::Any
	phi_eps::Any
	source::Any
end

# v, lambda and f for smooth steady state
function source_smooth(node)
	x = node[1]
	pi^2*cos(pi*x)
end

function v0_smooth(x)
	cos(pi*x)
end

function lambda0_smooth(x)
	if x > -1/2 && x < 1/2
		1.0
	else
		0.0
	end
end

# v, lambda and f for the rough steady state
function source_rough(node)
	x = node[1]
	if x > -1/2 && x < 1/2
		2*pi^2*cos(2*pi*x)
	else
		0.0
	end
end

function v0_rough(x)
	if x > -1/2 && x < 1/2
		1/2*cos(2*pi*x)
	else
		-1/2
	end
end


function lambda0_rough(x)
	if x > -1/4 && x < 1/4
		1.0
	else
		0.0
	end
end

function main(;type = "Smooth", visualize = true, test = false, p = 2)

	if type == "Smooth"
		source = source_smooth
	else
		source = source_rough
	end

	problem_data_PWL = ProblemData(
		ForwardBackwardHeat.PWLinear.phi, 
		ForwardBackwardHeat.PWLinear.kappa, 
		ForwardBackwardHeat.PWLinear.gamma_0, 
		ForwardBackwardHeat.PWLinear.phi_0, 
		source
	)
	H, lp_error_PWL, l2h1_error_PWL = convergence_error(problem_data_PWL; type)

	problem_data_Par = ProblemData(
		ForwardBackwardHeat.Parallel.phi, 
		ForwardBackwardHeat.Parallel.kappa, 
		ForwardBackwardHeat.Parallel.gamma_0, 
		ForwardBackwardHeat.Parallel.phi_0, 
		source
	)
	_, lp_error_Par, l2h1_error_Par = convergence_error(problem_data_Par; type)

	problem_data_Cub = ProblemData(
		ForwardBackwardHeat.Cubic.phi, 
		ForwardBackwardHeat.Cubic.kappa, 
		ForwardBackwardHeat.Cubic.gamma_0, 
		ForwardBackwardHeat.Cubic.phi_0, 
		source
	)
	_, lp_error_Cub, l2h1_error_Cub = convergence_error(problem_data_Cub; type)

	if visualize
		p = GridVisualizer(; Plotter = Plots, layout = (1, 2), legend = :rt, fast = true, size = (1000, 600))
		scalarplot!(
			p[1, 1], log10.(H.^(-1)), log10.(lp_error_PWL.^(-1)); label = "PWL", color = :blue,  
			title = "Convergence error in L^2_{t,x}", xlabel = "-log(dx)", ylabel = "-log-Error",
			markershape = :circle, markevery = 1
		)
		scalarplot!(p[1, 1], log10.(H.^(-1)), log10.(lp_error_Par.^(-1)); label = "Par", color = :green,clear = false)
		scalarplot!(p[1, 1], log10.(H.^(-1)), log10.(lp_error_Cub.^(-1)); label = "Cub", color = :red, clear = false)
		scalarplot!(
			p[1, 2], log10.(H.^(-1)), log10.(l2h1_error_PWL.^(-1)); label = "PWL", color = :blue,
			title = "Convergence error in L^2_tH^1_x", xlabel = "-log(dx)", ylabel = "-log-Error", 
			markershape = :circle, markevery = 1
		)
		scalarplot!(p[1, 2], log10.(H.^(-1)), log10.(l2h1_error_Par.^(-1)); label = "Par", color = :green, clear = false)
		scalarplot!(p[1, 2], log10.(H.^(-1)), log10.(l2h1_error_Cub.^(-1)); label = "Cub", color = :red, clear = false)
		reveal(p)
	end
	
	if test
		print("ConvergenceError.jl successful")
	end

	H, [lp_error_PWL, lp_error_Par, lp_error_Cub], [l2h1_error_PWL, l2h1_error_Par, l2h1_error_Cub]
	
end

function convergence_error(pb_data; type = "Smooth", p = 2)

	if type == "Smooth"
		v0_ini = v0_smooth
		lambda0_ini = lambda0_smooth
	else
		v0_ini = v0_rough
		lambda0_ini = lambda0_rough
	end

	# Unpacking problem data
	gamma_0, kappa = pb_data.gamma_eps, pb_data.kappa

	# Space-time bounds
	X_min = -1
	X_max = 1

	x0 = (X_max+X_min)/2

	T_max = 1

	u0, v0, lambda0 = prepare_IC(
		gamma_0, kappa, 
		v0_ini, lambda0_ini
	)
	
	# Sequence of mesh sizes
	H =[0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
	lp_error = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	l2h1_error = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	for i = 1:length(H)
		h = H[i]

		# Space discretisation
		X = X_min:h:X_max
		sgrid = simplexgrid(X)

		# Time discretisation
		k = 0.0025
		T = k:k:T_max

		tu, tv, tlambda, sys = fbheat(
			u0, v0, lambda0, 
			sgrid, T, pb_data
		)

		lp_er = 0
		l2h1_er = 0

		for j in 1:(length(tu.t)-1)
			tj = tu.t[j]
			tjp = tu.t[j + 1]
			dt = tjp - tj
			lp_er = lp_er + dt*(
				VoronoiFVM.lpnorm(sys, tv(tj) - tv(k), p)^p
				+ VoronoiFVM.lpnorm(sys, tu(tj) - tu(k), p)^p
				+ VoronoiFVM.lpnorm(sys, tlambda(tj) - tlambda(k), p)^p
			)
			l2h1_er = l2h1_er + dt*VoronoiFVM.h1norm(sys, tv(tj) - tv(k))^2
		end

		lp_error[i] = (lp_er)^(1/p)
		l2h1_error[i] = (l2h1_er)^(1/2)
	end

	H, lp_error, l2h1_error

end

end