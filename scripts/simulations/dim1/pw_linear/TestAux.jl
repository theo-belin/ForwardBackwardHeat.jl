module TestAux

using Printf
using PyPlot
using ForwardBackwardHeat
using ForwardBackwardHeat.PWLinear
using GridVisualize

function main()
    h = 0.1
    U = -7:h:7
    V = -5:h:5

    p = GridVisualizer(;Plotter = PyPlot, layout = (1, 4))
    scalarplot!(
        p[1,1], U, phi.(U);
        color = :blue,
        markershape = :circle, markevery = 1
    )
    scalarplot!(
        p[1, 2], U, phi_0.(U, 1);
        color = :blue,
        markershape = :circle, markevery = 1
    )
    scalarplot!(
        p[1, 3], U, phi_0.(U, 0.5);
        color = :blue,
        markershape = :circle, markevery = 1
    )
    scalarplot!(
        p[1, 4], U, phi_0.(U, 0);
        color = :blue,
        markershape = :circle, markevery = 1
    )
    reveal(p)
end

end
