# ForwardBackwardHeat

This project is based on the julia package [VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl) for finite volume discretisation of reaction-advection-diffusion equations.

[![Build Status](https://github.com/j-fu/ForwardBackwardHeat.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/j-fu/ForwardBackwardHeat.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Installing Julia

To install Julia, follow the Julia documentation at https://julialang.org/downloads/.

# Installing libraries and activating the project

Open your terminal, go to the local folder containing the files of the project and run 
```
% julia --project=.
```
and then in the REPL run
```
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.resolve()
```
This installs the julia packages needed to run the code. In case the general registry of packages available in Julia is not yet installed, then you should first run
```
julia> Pkg.Registry.add("General")
```

# Running a simulation

Choose a script, e.g. 'simulations/dim1/pw_linear/MovingFronts.jl' and run
```
julia> include("simulations/dim1/pw_linear/MovingFronts.jl")
julia> MovingFronts.main()
```
(The first time, running 'Module.main()' can take a while (depending on your processing power), this is due to the precompilation step embedded in Julia. The following runs are usually much faster.)

# Modify a script or files in src

You can modify a script or the source code directly, just save the modifications. They should take effect in the next run, without the need to reload the project, thanks to the package Revise.jl.
