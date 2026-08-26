c0 = -1
c1 = 1
P0 = -1
P1 = 1

k0 = 1
k1 = 1

function phi0(u)
    k0*(u - c0) + P1
end

function phi1(u)
    k1*(u - c1) + P0
end

function phi(u)
    if u <= c0
        phi0(u)
    elseif u > c0 && u < c1
        - ((P1 - P0)/(c1 - c0))*u
    elseif u >= c1
        phi1(u)
    end
end

function kappa0(v)
    (v - P1)/k0 + c0
end

function kappa1(v)
    (v - P0)/k1 + c1
end

kappa(v, lambda) = (1-lambda)*kappa0(v) + lambda*kappa1(v)
gamma_0 = construct_gamma_eps(P0, P1; eps = 0)
gamma_eps = construct_gamma_eps(P0, P1; eps = 0.01)

s0 = kappa0(P0)
s1 = kappa1(P1)

function phi_0(u,lambda)
    if u <= s0
        phi(u)
    elseif u > s0 && u < s0 + (c1 - s0)*lambda
        P0
    elseif u >= s0 + (c1 - s0)*lambda && u <= s1 + (c0 - s1)*(1-lambda)
        k_inter = (P1 - P0)/(s1 + (c0 - s1)*(1-lambda) - (s0 + (c1 - s0)*lambda))
        k_inter*(u - (s0 + (c1 - s0)*lambda)) + P0
    elseif u > s1 + (c0 - s1)*(1-lambda) && u < s1
        P1
    else
        phi(u)
    end
end

export phi, c0, c1, P0, P1, kappa, gamma_0, phi_0