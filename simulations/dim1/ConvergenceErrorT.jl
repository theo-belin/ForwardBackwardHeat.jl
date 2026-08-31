module ConvergenceErrorT

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
function u_exact(kappa,x)
	kappa0(v) = kappa(v,0)
	kappa1(v) = kappa(v,1)
	if x > -1/2 && x < 1/2
		(kappa0(cos(pi*x)) + kappa1(0))/2
	else
		(kappa0(0) + kappa1(cos(pi*x)))/2
	end
end

# exact lambda solution
function lambda_exact(kappa,t,x)
	kappa0(v) = kappa(v,0)
	kappa1(v) = kappa(v,1)
	(u_exact(kappa, x) - kappa0(v_exact(t,x)))/(kappa1(v_exact(t,x)) - kappa0(v_exact(t,x)))
end

function main(; visualize = true, test = false, p = 2)

	problem_data_PWL = ProblemData(
		ForwardBackwardHeat.PWLinear.phi, 
		ForwardBackwardHeat.PWLinear.kappa, 
		ForwardBackwardHeat.PWLinear.gamma_0, 
		ForwardBackwardHeat.PWLinear.phi_0, 
		source
	)
	K, lp_error_v_PWL, lp_error_u_PWL, l2h1_error_PWL = convergence_error(problem_data_PWL)

	problem_data_Par = ProblemData(
		ForwardBackwardHeat.Parallel.phi, 
		ForwardBackwardHeat.Parallel.kappa, 
		ForwardBackwardHeat.Parallel.gamma_0, 
		ForwardBackwardHeat.Parallel.phi_0, 
		source
	)
	_, lp_error_v_Par, lp_error_u_Par, l2h1_error_Par = convergence_error(problem_data_Par)

	problem_data_Cub = ProblemData(
		ForwardBackwardHeat.Cubic.phi, 
		ForwardBackwardHeat.Cubic.kappa, 
		ForwardBackwardHeat.Cubic.gamma_0, 
		ForwardBackwardHeat.Cubic.phi_0, 
		source
	)
	_, lp_error_v_Cub, lp_error_u_Cub, l2h1_error_Cub = convergence_error(problem_data_Cub)

	if visualize
		p = GridVisualizer(; Plotter = Plots, layout = (1, 3), legend = :rt, fast = true, size = (1000, 600))

		scalarplot!(
			p[1, 1], log10.(K.^(-1)), log10.(lp_error_v_PWL.^(-1)); label = "PWL", color = :blue, 
			title = "Convergence error of v in L^2_{t,x}", xlabel = "-log(dt)", ylabel = "-log-Error",
			markershape = :circle, markevery = 1
		)
		scalarplot!(p[1, 1], log10.(K.^(-1)), log10.(lp_error_v_Par.^(-1)); label = "Par", color = :green,clear = false)
		scalarplot!(p[1, 1], log10.(K.^(-1)), log10.(lp_error_v_Cub.^(-1)); label = "Cub", color = :red, clear = false)

		scalarplot!(
			p[1, 2], log10.(K.^(-1)), log10.(lp_error_u_PWL.^(-1)); label = "PWL", color = :blue, 
			title = "Convergence error of u in L^2_{t,x}", xlabel = "-log(dt)", ylabel = "-log-Error",
			markershape = :circle, markevery = 1
		)
		scalarplot!(p[1, 2], log10.(K.^(-1)), log10.(lp_error_u_Par.^(-1)); label = "Par", color = :green,clear = false)
		scalarplot!(p[1, 2], log10.(K.^(-1)), log10.(lp_error_u_Cub.^(-1)); label = "Cub", color = :red, clear = false)

		scalarplot!(
			p[1, 3], log10.(K.^(-1)), log10.(l2h1_error_PWL.^(-1)); label = "PWL", color = :blue,
			title = "Convergence error of v in L^2_tH^1_x", xlabel = "-log(dt)", ylabel = "-log-Error", 
			markershape = :circle, markevery = 1
		)
		scalarplot!(p[1, 3], log10.(K.^(-1)), log10.(l2h1_error_Par.^(-1)); label = "Par", color = :green, clear = false)
		scalarplot!(p[1, 3], log10.(K.^(-1)), log10.(l2h1_error_Cub.^(-1)); label = "Cub", color = :red, clear = false)

		reveal(p)
	end
	
	if test
		print("ConvergenceErrorT.jl successful")
	end

	K, [lp_error_PWL, lp_error_Par, lp_error_Cub], [l2h1_error_PWL, l2h1_error_Par, l2h1_error_Cub]
	
end

function convergence_error(pb_data; p = 2)

	# Unpacking problem data
	gamma_0, kappa = pb_data.gamma_eps, pb_data.kappa

	# Space-time bounds
	X_min = -1
	X_max = 1

	x0 = (X_max+X_min)/2

	T_max = 1

	v0(x) = v_exact(0,x)
	u0(x) = u_exact(kappa,x)
	lambda0(x) = lambda_exact(kappa,0,x)
	
	# Sequence of time steps
	K =[0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
	lp_error_v = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	lp_error_u = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	l2h1_error = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	for i = 1:length(K)
		k = K[i]

		# Space discretisation
		X = X_min:0.002:X_max
		sgrid = simplexgrid(X)

		# Time discretisation
		T = k:k:T_max

		tu, tv, tlambda, sys = fbheat(
			u0, v0, lambda0, 
			sgrid, T, pb_data
		)

		lp_er_v = 0
		lp_er_u = 0
		l2h1_er = 0

		tu_ex = unknowns(sys)
		tv_ex = unknowns(sys)
		tlambda_ex = unknowns(sys)
		tu_ex[1,:] = map(u0, sgrid)

		for j = 1:(length(tu.t)-1)
			tj = tu.t[j]
			tjp = tu.t[j + 1]
			dt = tjp - tj

			tv_ex[1,:] = map(x -> v_exact(tj, x), sgrid)
			tlambda_ex[1,:] = map(x -> lambda_exact(kappa, tj, x), sgrid)

			lp_er_v = lp_er_v + dt*(
				VoronoiFVM.lpnorm(sys, tv(tj) - tv_ex, p)^p
			)
			lp_er_u = lp_er_u + dt*(
				VoronoiFVM.lpnorm(sys, tu(tj) - tu_ex, p)^p
			)
				#+ VoronoiFVM.lpnorm(sys, tlambda(tj) - tlambda_ex, p)^p
			
			l2h1_er = l2h1_er + dt*VoronoiFVM.h1norm(sys, tv(tj) - tv_ex)^2
		end

		lp_error_v[i] = (lp_er_v)^(1/p)
		lp_error_u[i] = (lp_er_u)^(1/p)
		l2h1_error[i] = (l2h1_er)^(1/2)
	end

	K, lp_error_v, lp_error_u, l2h1_error

end

end