
function construct_gamma_eps(
    P0, P1; 
    eps = 0.1
    )
    """
        P0, P1: two real values.
            Required: P_1 > P_0.
        eps: small nonnegative value.

            Construct the (single-valued) function approximating the hysteresis operator. 
    """

    if eps > 0

        function gamma_eps(v, lambda)
            if v <= P0 - eps
                0
            elseif P0 - eps < v && v < P0 + eps*(2*lambda-1)
                (v - (P0 - eps))/(2*eps)
            elseif P0 + eps*(2*lambda-1) <= v && v <= P1 + eps*(2*lambda-1)
                lambda
            elseif P1 + eps*(2*lambda-1) < v && v < P1 + eps
                (v-(P1-eps))/(2*eps)
            else
                1
            end
        end
        gamma_eps

    else

        function gamma_0(v, lambda)
            if v < P0
                0
            elseif P0 <= v && v <= P1
                lambda
            else
                1
            end
        end
        gamma_0
    end
end

function construct_pwl(
    c0, c1, P0, P1; 
    k = 1, eps = 0.1
    )
    """
        c0, c1: local maximum and minimum (resp.) of the piecewise linear phi.
            Required: c0 < c1.
        P0, P1: local minimal and maximal value of phi.
            Required: P0 < P1.
        k: positive slope of phi below c0 and above c1.
        eps: regularization parameter, small, nonnegative.
        
            Computes: 
                - the unique piecewise linear function phi interpolating 
                    phi(c0) = P1 and phi(c1) = P0
                    with slope k below c0 and above c1,
                - the inverses kappa0, kappa1 of the increasing branches of phi,
                - the regularization phi_eps as the inverse of 
                    v |-> kappa(v, gamma_eps(v, lambda)).
    """
    function phi(u)
        if u <= c0
            k*(u-c0) + P1
        elseif c0 < u < c1
            (P1-P0)/(c0 - c1)*(u - c0) + P1
        else c1 <= u
            k*(u-c1) + P0
        end
    end

    kappa0(v) = (v-P1)/k + c0
    kappa1(v) = (v-P0)/k + c1

    delta = c1 - c0 + (P1 - P0)/k

    gamma_eps = construct_gamma_eps(P0, P1; eps)
    
    function phi_eps(u, lambda)
        if u <= c0 - (P1 - P0)/k - eps/k
            phi(u)
        elseif c0 - (P1 - P0)/k - eps/k < u < c0 - (P1 - P0)/k + (2*lambda -1)*eps/k + lambda*delta
        elseif c0 - (P1 - P0)/k + (2*lambda -1)*eps/k + lambda*delta <= u <= c1 + (P1 - P0)/k + (2*lambda -1)*eps/k - (1-lambda)*delta
            k*(u - c0 - lambda*delta) + P1
        elseif c1 + (P1 - P0)/k + (2*lambda -1)*eps/k - lambda*delta < u < c1 + (P1 - P0)/k - eps/k

        else
            phi(u)
        end
    end
    phi, kappa0, kappa1, phi_eps
end

function linrange_inserts(inserts, step)
    """ 
        inserts: a list of float values x0 < x1 < ... < xn.
            Required: n >= 1.
        step: a float.

            Computes a list of float values: [
                x0, x0 + step, ... x0 + k0*step, 
                x1, x1 + step, ..., x1 + k1*step, 
                ..., 
                x_(n-1),..., x_(n-1) + k_(n-1)*step
                xn].
    """

    n =  length(inserts)
    if n <= 1
        error("Error, there must be at least two elements in the list of inserts")
    else
        r = []
        for i = 1:n-1
            r = vcat(r, collect(inserts[i]:step:inserts[i + 1]))
        end
        if last(inserts) > last(r)
            push!(r, last(inserts))
        end
    end
    r
end

function invert_monotone(f, x = []; 
    xmin = 0, xmax = 1, special_vals = [], 
    step = 0.01, fill = 1)
    """
        f: a real monotone function to invert
        x: a sample of real arguments
        xmin, xmax: minimal value and maximal value of the domain of f
        special_vals: a list of particular point that have to be present in the sample
        step: the spacing of the sample
        fill: constant value for the extrapolation outside of the domain

            Computes the inverse of a real monotone function.
    """
    if x == []
        inserts = vcat(xmin, special_vals, xmax)
        x_vec = linrange_inserts(inserts, step)
    else
        x_vec = x[:]
    end

    y = f.(x_vec)
    Interpolations.deduplicate_knots!(y; move_knots = true)
    inv = Interpolations.linear_interpolation(y,x_vec, extrapolation_bc = fill)
    inv
end

function construct_kappa(phi, c0, c1; 
    umin = -5, umax = 5, special_vals = [], step = 0.01)
    """
        phi(u): real function
            - increasing from - oo to c0
            - decreasing from c0 to c1
            - increasing from c1 to +oo
        c0 < c1: real values
        umin, umax: range of values for the inverting
        special_vals: special values of interest (e.g. near critical points of phi)
        step: the discretization step for the inverting

            Constructs kappa0 and kappa1, the inverse of the two increasing branches of phi.
    """
    kappa0 = invert_monotone(
        phi; 
        xmin = umin, xmax = c0, 
        special_vals, 
        step, 
        fill = c0
        )

    kappa1 = invert_monotone(
        phi; 
        xmin = c1, xmax = umax, 
        special_vals, 
        step, 
        fill = c1)

    kappa(v,lambda) = (1-lambda)*kappa0(v) + lambda*kappa1(v)

    kappa
end

function construct_phi_eps(gamma_eps, kappa; 
    vmin = -5, vmax = 5, step_v = 0.01, step_lambda = 0.01)
    """
        kappa: (v, lambda) real valued function, monotone with respect to v and lambda
        gamma_eps: (v, lambda) real valued function
            Required:
                - between 0 and 1.
                - monotone with respect to v and lambda.
        v_min, v_max: range of the v values.
        step_v: discretization step for the v values for inverting.
        step_lambda: discretization step for the lambda values.

            Constructs phi_eps by inverting v |-> kappa(v,gamma_eps(v,lambda)) for a discrete set of lambdas
    """
    
    v_range = vmin:step_v:vmax
    phi_eps_list = [invert_monotone(
        v -> kappa(v,gamma_eps(v,lambda)), v_range;
        step = step_v
        ) for lambda in 0:step_lambda:1]

    function phi_eps(u,lambda)
        lambda_proj = max(min(1, lambda), 0)
        i = Int(round(lambda_proj/step_lambda)) + 1
        phi_eps_list[i](u)
    end
    phi_eps
end

function auxiliaries(phi, c0, c1, P0, P1;
    umin = -5, umax = 5,
    vmin = -5, vmax = 5,
    step = 0.01, eps = 0.01)
    # construction of the nondecreasing inverse branches
    kappa = construct_kappa(
        phi, c0, c1; 
        umin, umax,
        step)

    # Gamma_0 and Phi_0 (eps = 0)
    gamma_0 = construct_gamma_eps(P0, P1; eps = 0)
    phi_0 = construct_phi_eps(
        gamma_0, kappa;
        vmin, vmax)

    # Gamma_eps and Phi_eps (eps = 0.01)
    step_v = eps/10
    gamma_eps = construct_gamma_eps(P0, P1; eps)
    phi_eps = construct_phi_eps(
        gamma_eps, kappa; 
        vmin, vmax,
        step_v)

    kappa, gamma_0, phi_0, gamma_eps, phi_eps
end
