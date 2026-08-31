module ConvergenceErrorTX

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

# time-dependent coefficient
function A(t)
	m = 0.2 #should be chosen small enough to ensure time independence of u
	m*t^2
end

# time-dependent source term f
function source(node)
	x = node[1]
	t = node.time
	A(t)*pi^2*cos(pi*x)
end

# exact v solution
function v_exact(t,x)
	A(t)*cos(pi*x)
end

#exact u solutions
function u_exact(kappa,t,x)
	kappa0(v) = kappa(v,0)
	kappa1(v) = kappa(v,1)
	if x > -1/2 && x < 1/2
		#(kappa0(cos(pi*x)) + kappa1(0))/2
		kappa1(0)
	else
		#(kappa0(0) + kappa1(cos(pi*x)))/2
		kappa0(0)
	end
end

# exact lambda solution
function lambda_exact(kappa,t,x)
	kappa0(v) = kappa(v,0)
	kappa1(v) = kappa(v,1)
	(u_exact(kappa, t, x) - kappa0(v_exact(t, x)))/(kappa1(v_exact(t, x)) - kappa0(v_exact(t, x)))
end

function main(; visualize = true, test = false, p = 2)

	problem_data_PWL = ProblemData(
		ForwardBackwardHeat.PWLinear.phi, 
		ForwardBackwardHeat.PWLinear.kappa, 
		ForwardBackwardHeat.PWLinear.gamma_0, 
		ForwardBackwardHeat.PWLinear.phi_0, 
		source
	)
	K, _, lp_t_v_PWL, lp_t_u_PWL, lp_t_lambda_PWL, l2h1_t_PWL = convergence_error(problem_data_PWL; type = "time")
	_, H, lp_x_v_PWL, lp_x_u_PWL, lp_x_lambda_PWL, l2h1_x_PWL = convergence_error(problem_data_PWL; type = "space")

	problem_data_Par = ProblemData(
		ForwardBackwardHeat.Parallel.phi, 
		ForwardBackwardHeat.Parallel.kappa, 
		ForwardBackwardHeat.Parallel.gamma_0, 
		ForwardBackwardHeat.Parallel.phi_0, 
		source
	)
	_, _, lp_t_v_Par, lp_t_u_Par, lp_t_lambda_Par, l2h1_t_Par = convergence_error(problem_data_Par; type = "time")
	_, _, lp_x_v_Par, lp_x_u_Par, lp_x_lambda_Par, l2h1_x_Par = convergence_error(problem_data_Par; type = "space")

	problem_data_Cub = ProblemData(
		ForwardBackwardHeat.Cubic.phi, 
		ForwardBackwardHeat.Cubic.kappa, 
		ForwardBackwardHeat.Cubic.gamma_0, 
		ForwardBackwardHeat.Cubic.phi_0, 
		source
	)
	_, _, lp_t_v_Cub, lp_t_u_Cub, lp_t_lambda_Cub, l2h1_t_Cub = convergence_error(problem_data_Cub; type = "time")
	_, _, lp_x_v_Cub, lp_x_u_Cub, lp_x_lambda_Cub, l2h1_x_Cub = convergence_error(problem_data_Cub; type = "space")

	if visualize
		p = GridVisualizer(; Plotter = Plots, layout = (2, 4), legend = :rt, fast = true, size = (1000, 600))

		labels = ["PWL", "Par", "Cub"]
		colors = [:blue, :green, :red]

		plot_diagram(p[1, 1], K, [lp_t_v_PWL, lp_t_v_Par, lp_t_v_Cub], "Convergence error of v in L^2_{t,x}, with dt", "-log(dt)", "-log-Error", labels, colors)
		plot_diagram(p[1, 2], K, [lp_t_u_PWL, lp_t_u_Par, lp_t_u_Cub], "Convergence error of u in L^2_{t,x}, with dt", "-log(dt)", "-log-Error", labels, colors)
		plot_diagram(p[1, 3], K, [lp_t_lambda_PWL, lp_t_lambda_Par, lp_t_lambda_Cub], "Convergence error of lambda in L^2_{t,x}, with dt", "-log(dt)", "-log-Error", labels, colors)
		plot_diagram(p[1, 4], K, [l2h1_t_PWL, l2h1_t_Par, l2h1_t_Cub], "Convergence error of v in L^2_tH^1_x, with dt", "-log(dt)", "-log-Error", labels, colors)
		plot_diagram(p[2, 1], H, [lp_x_v_PWL, lp_x_v_Par, lp_x_v_Cub], "Convergence error of v in L^2_{t,x}, with dx", "-log(dx)", "-log-Error", labels, colors)
		plot_diagram(p[2, 2], H, [lp_x_u_PWL, lp_x_u_Par, lp_x_u_Cub], "Convergence error of u in L^2_{t,x}, with dx", "-log(dx)", "-log-Error", labels, colors)
		plot_diagram(p[2, 3], H, [lp_x_lambda_PWL, lp_x_lambda_Par, lp_x_lambda_Cub], "Convergence error of lambda in L^2_{t,x}, with dx", "-log(dx)", "-log-Error", labels, colors)
		plot_diagram(p[2, 4], H, [l2h1_x_PWL, l2h1_x_Par, l2h1_x_Cub], "Convergence error of v in L^2_tH^1_x, with dx", "-log(dx)", "-log-Error", labels, colors)
		reveal(p)
	end
	
	if test
		print("ConvergenceErrorT.jl successful")
	end

	K, [lp_t_v_PWL, lp_t_u_PWL, lp_t_lambda_PWL, l2h1_t_PWL], [lp_t_v_Par, lp_t_u_Par, lp_t_lambda_Par, l2h1_t_Par], [lp_t_v_Cub, lp_t_u_Cub, lp_t_lambda_Cub, l2h1_t_Cub], H, [lp_x_v_PWL, lp_x_u_PWL, lp_x_lambda_PWL, l2h1_x_PWL], [lp_x_v_Par, lp_x_u_Par, lp_x_lambda_Par, l2h1_x_Par], [lp_x_v_Cub, lp_x_u_Cub, lp_x_lambda_Cub, l2h1_x_Cub]
	
end

function plot_diagram(p, steps, errors, title_str, xlabel_str, ylabel_str, labels, colors)
	scalarplot!(
			p, log10.(steps.^(-1)), log10.(errors[1].^(-1)); title = title_str, xlabel = xlabel_str, ylabel = ylabel_str,
			label = labels[1], color = colors[1], markershape = :circle, markevery = 1
		)
	for i = 2:length(errors)
		scalarplot!(
			p, log10.(steps.^(-1)), log10.(errors[i].^(-1)); label = labels[i], color = colors[i], markershape = :circle, markevery = 1
		)
	end
end

function compute_errors(v_ex, u_ex, lambda_ex, sgrid, T, pb_data; p = 2)

	lp_error_v = 0
	lp_error_u = 0
	lp_error_lambda = 0
	l2h1_error = 0

	v0(x) = v_ex(0,x)
	u0(x) = u_ex(0,x)
	lambda0(x) = lambda_ex(0,x)

	tu, tv, tlambda, sys = fbheat(
		u0, v0, lambda0, 
		sgrid, T, pb_data
	)

	tu_ex = unknowns(sys)
	tv_ex = unknowns(sys)
	tlambda_ex = unknowns(sys)

	for j = 1:(length(tu.t)-1)
		tj = tu.t[j]
		tjp = tu.t[j + 1]
		dt = tjp - tj

		tv_ex[1,:] = map(x -> v_ex(tj, x), sgrid)
		tu_ex[1,:] = map(x -> u_ex(tj, x), sgrid)
		tlambda_ex[1,:] = map(x -> lambda_ex(tj, x), sgrid)

		lp_error_v += dt * (VoronoiFVM.lpnorm(sys, tv(tj) - tv_ex, p)^p)
		lp_error_u += dt * (VoronoiFVM.lpnorm(sys, tu(tj) - tu_ex, p)^p)
		lp_error_lambda += dt * (VoronoiFVM.lpnorm(sys, tlambda(tj) - tlambda_ex, p)^p)
		l2h1_error += dt * VoronoiFVM.h1norm(sys, tv(tj) - tv_ex)^2
	end

	lp_error_v^(1/p), lp_error_u^(1/p), lp_error_lambda^(1/p), l2h1_error^(1/2)
end

function convergence_error(pb_data; type = "time", p = 2)

	# Unpacking problem data
	gamma_0, kappa = pb_data.gamma_eps, pb_data.kappa


	# Space-time bounds
	X_min = -1
	X_max = 1

	x0 = (X_max+X_min)/2

	T_max = 1
	
	if type == "time"
		# Sequence of time steps
		K =[0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
		H = [0.001] # fixed mesh size for time convergence
	else
		# Sequence of mesh sizes
		K = [0.001] # fixed time step for space convergence
		H = [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
	end
	
	lp_error_v = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	lp_error_u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	lp_error_lambda = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	l2h1_error = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	u_ex(t,x) = u_exact(kappa,t,x)
	lambda_ex(t,x) = lambda_exact(kappa,t,x)

	i = 1
	for k = K
		for h = H
			# Space discretisation
			X = X_min:h:X_max
			sgrid = simplexgrid(X)
			# Time discretisation
			T = k:k:T_max

			lp_er_v, lp_er_u, lp_er_lambda, l2h1_er = compute_errors(v_exact, u_ex, lambda_ex, sgrid, T, pb_data; p = p)

			lp_error_v[i] = lp_er_v
			lp_error_u[i] = lp_er_u
			lp_error_lambda[i] = lp_er_lambda
			l2h1_error[i] = l2h1_er
			i += 1
		end
	end

	K, H, lp_error_v, lp_error_u, lp_error_lambda, l2h1_error

end

end