ref(u) = sin(pi*u)/(2*pi) + u
c0 = -1.0
c1 = 1.0
P0 = -1.0
P1 = 1.0

function phi(u)
    if u <= c0
        ref(u + 2)
    elseif u > c0 && u < c1
        - u
    elseif u >= c1
        ref(u - 2)
    end
end

gamma_0 = construct_gamma_eps(P0, P1; eps = 0)
gamma_eps = construct_gamma_eps(P0, P1; eps = 0.01)

inv_ref = invert_monotone(ref; xmin = -5, xmax = 5, step = 0.001)

function kappa0(v)
    inv_ref(v) - 2
end

function kappa1(v)
    inv_ref(v) + 2
end

function kappa(v,lambda)
    (1-lambda)*kappa0(v) + lambda*kappa1(v)
end

s0 = kappa0(P0)
s1 = kappa1(P1)

function phi_0(u,lambda)
    if u <= s0
        phi(u)
    elseif u > s0 && u < s0 + 4*lambda
        P0
    elseif u >= s0 + 4*lambda && u <= s1 - 4*(1-lambda)
        ref(u + 4*(1/2 - lambda))
    elseif u > s1 - 4*(1-lambda) && u < s1
        P1
    else
        phi(u)
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

export phi, c0, c1, P0, P1, kappa, gamma_0, phi_0