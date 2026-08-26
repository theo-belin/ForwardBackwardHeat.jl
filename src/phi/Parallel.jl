ref(u) = sin(pi*u)/(2*pi) + u
c0 = -1.0
c1 = 1.0
P0 = -1.0
P1 = 1.0

function phi(u)
    if u <= c0
        ref(u + 2)
    elseif u > c0 && u <= c1
        - u
    elseif u >= c1
        ref(u - 2)
    end
end

# function phi_0(u,lambda)
#     if u <= s0
#         ref(u + 1)
#     elseif u > s0 && u < s0 + 2*lambda
#         P0
#     elseif u >= s0 + 2*lambda && u <= s1 - 2*(1-lambda)
#         ref(u + 1 - 2*lambda)
#     elseif u > s1 - 2*(1-lambda) && u < s1
#         P1
#     else
#         ref(u - 1)
#     end
# end

kappa, gamma_0, phi_0, gamma_eps, phi_eps = auxiliaries(
    phi, c0, c1, P0, P1;
    umin = -7, umax = 7,
    vmin = -5, vmax = 5,
    step = 0.001)

export phi, c0, c1, P0, P1, kappa, gamma_0, phi_0, gamma_eps, phi_eps