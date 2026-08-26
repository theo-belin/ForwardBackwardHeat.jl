module TestAux

using Printf
using PyPlot
using ForwardBackwardHeatProject
using ForwardBackwardHeatProject.Cubic
using GridVisualize

function main()
    h = 0.1
    U = -1.5:h:1.5
    V = (P0-2):h:(P1 + 2)
    L = [0, 0.1, 0.5, 1]

    for i = (1:length(L))
        p = GridVisualizer(;Plotter = PyPlot, layout = (2, 3))
        scalarplot!(
            p[1, 1], U, phi_0.(U, L[i]); title = L"$\Phi_0$",
            color = :blue, label = "numerical",
            markershape = :circle, markevery = 1
        )
        scalarplot!(
            p[1, 2], U, phi_eps.(U, L[i]); title = L"$\Phi_\epsilon$",
        color = :blue, label = "numerical",
        markershape = :circle, markevery = 1
        )
        scalarplot!(
            p[2, 1], V, kappa.(V, L[i]); title = L"$\kappa$",
            color = :blue, label = "numerical",
            markershape = :circle, markevery = 1
        )
        scalarplot!(
            p[2, 2], V, gamma_0.(V, L[i]); title = L"$\Gamma_0$",
            color = :blue, label = "numerical",
            markershape = :circle, markevery = 1
        )
        scalarplot!(
            p[2, 3], V, gamma_eps.(V, L[i]); title = L"$\Gamma_\epsilon$",
            color = :blue, label = "numerical",
            markershape = :circle, markevery = 1
        )
        reveal(p)
    end
end

end
