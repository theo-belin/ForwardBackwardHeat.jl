
function exact_riemann(ul, ur, ml, mr, c0, c1, P0, P1)
    """
        u_l
        u_r
        k_l
        k_r
        c0
        c1
        P0
        P1
    """
    
    #Exact solution on the whole line, without boundary conditions
    kl = P1 - ml*c0
    kr = P0 - mr*c1

    function phi(u)
        if u <= c0
            ml*u + kl
        elseif u >= c1
            mr*u + kr
        else
            error("The data is in the instable region !")
        end
    end

    function e_m_plus(m, xi)
        1/2*erfc(xi/sqrt(4*m))
    end

    function e_m_minus(m, xi)
        1/2*(1 + erf(xi/sqrt(4*m)))
    end

    function g(xi)
        fac = (sqrt(mr) - sqrt(ml))/(sqrt(mr) + sqrt(ml))
        jump_phi = phi(ur) - phi(ul)
        if xi < 0
            phi(ul)*e_m_plus(ml, xi) + phi(ur)*e_m_minus(ml, xi) - fac*jump_phi*e_m_minus(ml, xi)
        else
            phi(ul)*e_m_plus(mr, xi) + phi(ur)*e_m_minus(mr, xi) - fac*jump_phi*e_m_minus(mr, xi)
        end
    end

    function f(xi)
        
        if xi < 0
            (g(xi) - kl)/ml
        else
            (g(xi) - kr)/mr
        end
    end

    function u(t,x)
        f(x/sqrt(t))
    end
    u
end