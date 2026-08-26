function prepare_IC_2D(gamma, kappa, v0, lambda)
    """
        gamma: a function of v and lambda, which outputs the projection of lambda on the admissible graph of the hysteresis operator
        kappa: a function of v and lambda, which outputs the associated u value
        v0: a real function of x
        lambda: a real function of x
    """
    lambda0(x, y) = gamma(v0(x, y), lambda(x, y))
    u0(x, y) = kappa(v0(x, y), lambda0(x, y))

    u0, v0, lambda0
end

# Riemann-type initial conditions
function construct_2D_indicator_fun(shape; center = [0,0], size = 1)
    
    x0, y0 = center

    function indicator_fun(x,y)
        if shape == "HalfSpace"
            x < x0
        elseif shape == "Disk"
            (x-x0)^2 + (y-y0)^2 < size^2
        elseif shape == "Square"
            abs(x-x0) < size && abs(y-y0) < size
        else
            error("No indicator function matching given shape")
        end
    end

    indicator_fun
end


function IC_riemann_2D(gamma, kappa, 
    v_in, v_out, lambda_in, lambda_out, 
    indicator_fun)
    """
    gamma, kappa : parameters of the nonlinearity
    v_out, lambda_out : the values of (v, lambda) outside the shape 
    v_in, lambda_in : the values of (v, lambda) inside the shape
    indicator_fun : the indicator function of the shape (optional)
    """

    function v0(x,y)
        if indicator_fun(x,y)
            v_in
        else
            v_out
        end
    end

    function lambda(x,y)
        if indicator_fun(x,y)
            lambda_in
        else
            lambda_out
        end
    end

    u0, v0, lambda0 = prepare_IC_2D(gamma, kappa, v0, lambda)
end

## General initial conditions 

# Smooth Checker board
function IC_smooth_checker_board(
    gamma, kappa,
    amp, freq; lambda_prior = 1/2)
    """
    gamma, kappa : parameters of the nonlinearity
    amp, freq : amplitude and frequency of the checker board
    prior_lambda : the default value of the phase fraction when v0 is between P0 and P1
    """
    function v0(x,y)
        amp*sin(freq*x)*sin(freq*y)
    end

    lambda0(x,y) = prior_lambda

    u0, v0, lambda0 = prepare_IC_2D(gamma, kappa, v0, lambda0)
end

# Random initial condition
function IC_sampled_random_2D(gamma, kappa, x_min, x_max, y_min, y_max, gamma_eps; seed = nothing, v_min = -1, v_max = 1, samples = 100)
    """
    gamma_eps : to perform the projection of the phase fraction
    kappa : to define u0
    seed : setting the random seed for repeatability
    v_min, v_max : min, max values for v
    """

    if !isnothing(seed)
        Random.seed!(seed)
    end

    samp_v = v_min .+ (v_max-v_min).*rand(samples, samples)
    samp_lambda = rand(samples, samples)

    x_length = x_max - x_min
    y_length = y_max - y_min

    function v0(x, y)
        i = Int(round((x - x_min)*(samples-1)/x_length)) + 1
        j = Int(round((y - y_min)*(samples-1)/y_length)) + 1
        samp_v[i, j]
    end

    function lambda_rand(x, y)
        i = Int(round((x - x_min)*(samples-1)/x_length)) + 1
        j = Int(round((y - y_min)*(samples-1)/y_length)) + 1
        samp_lambda[i, j]
    end

    u0, v0, lambda0 = prepare_IC_2D(v0, lambda_rand)
end

function IC_mixed_2D(kappa, data_left, data_right; x0 = 0, y0 = 0, k = 1)
    """
    kappa : to define u0
    data_left : data on the left of the straight interface
    data_right : fixed data on the right of the straigh interface. A linear function is added to exhibit the two behaviours. 

    x0, y0: x-location of the initial interface, y location of the start of the linear perturbation.
    k: coefficient of the linear perturbation.
    """

    v_l, lambda_l = data_left
    v_r, lambda_r = data_right

    function v0(x,y)
        if x < x0
            v_l
        else
            if y < y0
                v_r
            else
                v_r + (y-y0)*k
            end
        end
    end

    function lambda0(x,y)
        if x < x0
            lambda_l
        else
            lambda_r
        end
    end

    u0(x,y) = kappa(v0(x,y), lambda0(x,y))

    u0, v0, lambda0
end