module L1Contraction

using Printf

using ForwardBackwardHeat
using ForwardBackwardHeat.Cubic
using ExtendableGrids, Plots
using GridVisualize
using VoronoiFVM

# Numerical Evidence for L^1 Translation Estimate

mutable struct ProblemData
    phi::Any
	kappa::Any
	gamma_eps::Any
	phi_eps::Any
	source::Any
end

function main(; visualize = true, test = false)

	# Spatial domain
	X_min = -1
	X_max = 1
	# Space discretization
	h = 0.2
	X=X_min:h:X_max
	sgrid = simplexgrid(X)

	# Time discretization
	T_max = 0.5
	k = h/40
	T=k:k:T_max

	# First test with saturated moving front
	x0 = (X_max+X_min)/2

	problem_data = ProblemData(phi, kappa, gamma_0, phi_0, _ -> 0)
	
	u0_mf, v0_mf, lambda0_mf = IC_riemann_1D(
		gamma_0, kappa,
		P1 + 4, P1, 
		1, 0;
		x0 = x0
	)

	tu_mf, tv_mf, tlambda_mf, sys = fbheat(
		u0_mf, v0_mf, lambda0_mf,  
		sgrid, T, problem_data
	)

	off_set = 0.1
	u0_mf_os, v0_mf_os, lambda0_mf_os = IC_riemann_1D(
		gamma_0, kappa, 
		P1 + 4, P1,
		1, 0; 
		x0 = x0 + off_set
		)

	tu_mf_os, tv_mf_os, tlambda_mf_os, _ = fbheat(
		u0_mf_os, v0_mf_os, lambda0_mf_os,
		sgrid, T, problem_data
	)

	# Construct the quantity that has L^1 contraction (implies L^1 bounds for lambda and v)
	contracting_fun = build_contracting_fun(
		tlambda_mf, tlambda_mf_os,
		tv_mf, tv_mf_os
	)

	# discrep = build_W11_Linf_discrepancy(
	# 	tv_mf, tv_mf_os, 
	# 	k, T[1:80:end]
	# )

	l1_contracting_fun = [VoronoiFVM.lpnorm(sys, contracting_fun(t), 1) for t = T]
	l1_u = [VoronoiFVM.lpnorm(sys, tu_mf(t) - tu_mf_os(t), 1) for t = T]
	#discrep = [integrated(t)*l1_norms[1] for t = T[1:40:end]]

	if visualize
		p = GridVisualizer(; Plotter = Plots, layout = (1, 2), fast = true)
		scalarplot!(p[1,1], T, l1_u .- l1_u[1]; title = "Evolution of the L^1 norm of u(t,cdot) - u(t, cdot + h)", size = (500, 200), xlabel = "t", ylabel = "|u(t) - u_h(t)|_{L^1} - |u_0 - u_{0, h}|_{L^1}")
		scalarplot!(p[1, 2], T, l1_contracting_fun .- l1_contracting_fun[1] ; title = "Evolution of the L^1 norm of the contracting function", size = (500, 200), xlabel = "t", ylabel = "|f(t)|_{L^1} - |f(0)|_{L^1}")
		reveal(p)
	end

	if test
		print("L1Contraction.jl successful")
	end
end

function build_contracting_fun(tlambda1, tlambda2, tv1, tv2)

	kappa0(v) = kappa(v, 0)
	kappa1(v) = kappa(v, 1)

	c_lambda(t) = (kappa1.(tv1(t)) .+  kappa1.(tv2(t)) .- (kappa0.(tv1(t))  .+ kappa0.(tv2(t))))./2
	c0_v(t) = (1 .- (tlambda1(t) .+ tlambda2(t))./2)
	c1_v(t) = 1 .- c0_v(t)
	fun(t) = c_lambda(t).*abs.(tlambda1(t) .- tlambda2(t)) + c0_v(t).*abs.(kappa0.(tv1(t)) - kappa0.(tv2(t))) .+ c1_v(t).*abs.(kappa1.(tv1(t)) .- kappa1.(tv2(t)))
	
	fun
end

function build_W11_Linf_discrepancy(tv1, tv2, k, T)
	c(t) = (kappa1.(tv1(t)) .+  kappa1.(tv2(t)) .- (kappa0.(tv1(t))  .+ kappa0.(tv2(t))))./2
	dtc_c(t) = maximum(abs.(c(min(t + k, T[end])) .- c(t))./c(t))

	vals = [dtc_c(t) for t = T]
	sums = zeros(size(vals, 1), 1)
	i = 2
	for v = vals[1:end-1]
		sums[i] = sums[i-1] + v*80
		i = i + 1
	end

	sums
end

end