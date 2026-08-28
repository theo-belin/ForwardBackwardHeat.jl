using Printf
using PyPlot
using ForwardBackwardHeat
using ForwardBackwardHeat.Parallel
using GridVisualize

function main()
    h = 0.1
    U = -7:h:7
    V = -5:h:5

    p = GridVisualizer(;Plotter = PyPlot, layout = (2, 4))
    scalarplot!(
        p[1,1], U, phi.(U);
        color = :blue
    )
    scalarplot!(
        p[1, 2], U, phi_0.(U, 1);
        color = :blue
    )
    scalarplot!(
        p[1, 3], U, phi_0.(U, 0.5);
        color = :blue
    )
    scalarplot!(
        p[1, 4], U, phi_0.(U, 0);
        color = :blue
    )
    scalarplot!(
        p[2,1], V, gamma_0.(V, 0);
        color = :blue
    )
    scalarplot!(
        p[2, 2], V, gamma_0.(V, 0.5);
        color = :blue
    )
    scalarplot!(
        p[2, 3], V, kappa.(V, 1);
        color = :blue
    )
    scalarplot!(
        p[2, 4], V, kappa.(V, 0);
        color = :blue
    )
    reveal(p)
end
