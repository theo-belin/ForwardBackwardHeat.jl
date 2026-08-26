# Riemann initial conditions
function prepare_IC(gamma, kappa, v0, lambda)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        v0: a real function of x
        lambda: a real function of x
    """
    lambda0(x) = gamma(v0(x), lambda(x))
    u0(x) = kappa(v0(x), lambda0(x))

    u0, v0, lambda0
end

function IC_riemann_1D(gamma, kappa, v_l, v_r, lambda_l, lambda_r; x0 = 0)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        v_l, v_r: value of v on the left and right of the discontinuity
        lambda_l, lambda_r: value of lambda on the left and right of the discontinuity
        x0: the position of the discontinuity 
    """

    function v0(x)
        if x < x0
            v_l
        else
            v_r
        end
    end

    function lambda(x)
        if x < x0
            lambda_l
        else
            lambda_r
        end
    end

    u0, v0, lambda0 = prepare_IC(gamma, kappa, v0, lambda)
end

# Sinusoidal initial condition
function IC_sinusoidal_1D(gamma, kappa; 
    amp = 1, 
    freq = 1, 
    off_set = 0, 
    lambda = x -> 1/2)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        amp, freq, off_set: amplitude, frequence and off_set of the sinusoidal signal
        lambda: a real function of x
    """
    function v0(x)
        amp*sin(freq*x) + off_set
    end

    u0, v0, lambda0 = prepare_IC(gamma, kappa, v0, lambda)
end

function IC_random_1D(gamma, kappa; 
    x_min = 0, x_max = 1, 
    v_min = -1, v_max = 1, 
    lambda = x -> 1/2, 
    sampling = 10, 
    seed = nothing)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        v_min, v_max: range of the random sampling for v
        sampling: number of points sampled
        lambda: a real function of x
    """
    if !isnothing(seed)
        Random.seed!(seed)
    end

    v_vec = v_min .+ (v_max-v_min).*rand(sampling)

    x_length = x_max - x_min

    function v0(x)
        i = Int(round((x - x_min)*(sampling-1)/x_length)) + 1
        v_vec[i]
    end

    u0, v0, lambda0 = prepare_IC(gamma, kappa, v0, lambda)
end

# Propagation front moving with a delay from left to right
function IC_pw_affine(gamma, kappa, 
    coeff_l_v, coeff_r_v, 
    lambda; 
    x0 = 0)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        coeff_l_v: pair of coefficients (a_l, b_l), where 
            a_l is the slope on the left
            b_l is the value at x0
        coeff_r_v: as above but for the right 
        lambda: a real function of x
        x0: the point separating the two slopes
    """
    a_l, b_l = coeff_l_v
    a_r, b_r = coeff_r_v

    function v0(x)
        if x < x0
            k_l*(x - x0) + b_l
        else
            k_r*(x - x0) + b_r
        end
    end

    u0, v0, lambda0 = prepare_IC(gamma, kappa, v0, lambda)
end