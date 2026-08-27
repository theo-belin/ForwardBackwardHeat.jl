module ForwardBackwardHeat

using VoronoiFVM, ExtendableGrids
using Interpolations
using SpecialFunctions
using Random
using Revise
import DrWatson

include("simulations.jl")

export fbheat

include("specialplots.jl")

export xtplot, videodir

include("initial_conditions/initial_conditions_1D.jl")

export IC_riemann_1D, IC_pw_affine, IC_sinusoidal_1D, IC_random_1D, prepare_IC

include("initial_conditions/initial_conditions_2D.jl")

export construct_2D_indicator_fun, IC_riemann_2D, u0_smooth_checker_board, IC_random_2D, IC_mixed_2D, IC_sampled_random_2D, prepare_IC_2D

include("function_build.jl")

export invert_monotone, construct_kappa_i, construct_gamma_eps, construct_phi_eps, auxiliaries

module Cubic
using ..ForwardBackwardHeat
include("phi/Cubic.jl")
end

module Parallel
using ..ForwardBackwardHeat
include("phi/Parallel.jl")
end

module PWLinear
using ..ForwardBackwardHeat
include("phi/PWLinear.jl")
end

end