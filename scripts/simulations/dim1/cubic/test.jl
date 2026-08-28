using Printf
using PyPlot
using ForwardBackwardHeat
using ExtendableGrids
using GridVisualize
using VoronoiFVM

include("L1Contraction.jl")
include("MovingFronts.jl")
include("MovingFrontsForcing.jl")
include("OscillatingIC.jl")
include("SteadyFront.jl")
include("SteadyState.jl")

function test()
    #1D tests
    L1Contraction.main(; visualize = false, test = true)
    MovingFronts.main(; visualize = false, test = true)
    MovingFrontsForcing.main(; visualize = false, test = true)
    OscillatingIC.main(; visualize = false, test = true)
    SteadyFront.main(; visualize = false, test = true)
    SteadyState.main(; visualize = false, test = true)
end

test()