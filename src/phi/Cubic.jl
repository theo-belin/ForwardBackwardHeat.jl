phi(u) = u^3 - u

c0 = -sqrt(3)/3
c1 = sqrt(3)/3
P1 = phi(c0)
P0 = phi(c1)

kappa, gamma_0, phi_0, gamma_eps, phi_eps = auxiliaries(phi, c0, c1, P0, P1)

export phi, c0, c1, P0, P1, kappa, gamma_0, phi_0, gamma_eps, phi_eps