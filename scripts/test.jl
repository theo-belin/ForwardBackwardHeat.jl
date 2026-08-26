using Printf
using PyPlot
using ForwardBackwardHeat
using ExtendableGrids
using GridVisualize
using VoronoiFVM

include("simulations/dim1/cubic/MovingFronts.jl")
include("simulations/dim1/cubic/MovingFrontsForcing.jl")
include("simulations/dim1/cubic/SteadyFront.jl")
include("simulations/dim1/cubic/OscillatingIC.jl")
include("simulations/dim1/cubic/ConvergenceError.jl")
include("simulations/dim1/cubic/L1Contraction.jl")
include("simulations/dim2/pw_linear/HalfSpaceFront.jl")
include("simulations/dim2/pw_linear/ExpandingDisk.jl")
include("simulations/dim2/pw_linear/TwoExpandingDisks.jl")

function test()
    MovingFronts.main(; visualize = false, test = true)
    MovingFrontsForcing.main(; visualize = false, test = true)
    SteadyFront.main(; visualize = false, test = true)
    OscillatingIC.main(; visualize = false, test = true)
    ConvergenceError.main(; visualize = false, test = true)
    L1Contraction.main(; visualize = false, test = true)
    HalfSpaceFront.main(;visualize = false, test = true)
    ExpandingDisk.main(;visualize = false, test = true)
    TwoExpandingDisk.main(;visualize = false, test = true)
end

test()