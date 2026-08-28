function xtplot(
    X,
    tsol;
    Plotter = GridVisualize.Plots,
    tscale = :identity,
    tlabel = "t",
    kwargs...,
)
    T = tsol.t
    if tscale == :log
        T = log10.(T)
        if tsol.t[1] ≈ 0.0
            T[1] = log10(tsol.t[2] / 2)
        end
        if tlabel == "t"
            tlabel = "log10(t)"
        end
    end

    (xmin, xmax) = extrema(X)
    (tmin, tmax) = extrema(T)
    xtaspect = 0.7 * (tmax - tmin) / (xmax - xmin)
    scalarplot(
        simplexgrid(X, T),
        vec(tsol);
        Plotter,
        aspect = 1 / xtaspect,
        ylabel = tlabel,
        kwargs...,
    )
end

#videodir(args...)=DrWatson.projectdir("videos", args...)

# function create_axis(title, 
# 	y_min, y_max, 
# 	y_vals, 
# 	y_tick_labels, y_label;
# 	x_min = 0, x_max = 1.1)
# 	"""
# 		Create stereotypical axis of rectangular horizontal shape for useful plots
# 	"""
# 	ytick_str = replace(string(y_vals), "[" => "", "]" => "")

# 	@pgf axis = PGFPlotsX.Axis(
# 		{
# 		height = raw"7cm",
# 		width = raw"14cm",
# 		title = title,
# 		tick_label_style=raw"{font=small}",
	
# 		xmin = x_min,
# 		xmax = x_max,
# 		xtick = raw"{0, 0.5, 1}",
# 		xticklabels = raw"{0,x_0 = frac{1}{2}, 1}",
# 		xlabel_style = raw"{at = {(xticklabel cs:1,0)}}",
	
# 		ymin = y_min,
# 		ymax = y_max,
# 		ytick = ytick_str,
# 		yticklabels = y_tick_labels,
	
# 		xlabel = raw"x",
# 		ylabel = y_label,
# 		ylabel_style = raw"{at = {(yticklabel cs:1,0)}, rotate = -90}",
	
# 		axis_x_line = "bottom",
# 		axis_y_line = "left",
# 		}
# 	)
# 	axis

# end

# function pgfplot_u_1D(
# 	tu, 
# 	u_min, u_max, 
# 	u_vals, 
# 	u_tick_labels, 
# 	plot_times, X, 
# 	title, file_name
# 	)

# 	axis = create_axis(title, 
# 	u_min, u_max, u_vals, u_tick_labels, 
# 	raw"u")

# 	#styles = ["loosely dotted", "loosely dashed", "loosely dashdotted", "loosely dashdotdotted", "solid"]
# 	marks = ["|", "square", "star", "o", "x"]
# 	colors = ["red", "red!80!blue", "red!60!blue", "red!40!blue", "red!20!blue", "blue"]
# 	legends = [raw"t = "*string(t) for t = plot_times]
# 	i = 1
	
# 	for t = plot_times
# 		p = tu(t + 0.0025)[1,:]
# 		push!(axis, @pgf PGFPlotsX.Plot(
# 			{
# 				color = colors[i],
# 				#"only marks",
# 				mark = marks[i],
# 				mark_repeat = 6,
# 				mark_options = raw"{scale = 1.3}",
# 				x = "x",
# 				y = "u",
# 			},
# 			Table([:x => X, :v => p])
# 			)
# 		)
# 		i = i + 1
# 	end

# 	push!(axis, @pgf PGFPlotsX.Legend(legends))
	
# 	pgfsave("graphs/pdf/"*file_name*".pdf", axis)
# 	pgfsave("graphs/tikz/"*file_name*".tikz", axis)
	
# end

# function pgfplot_v_1D(
# 	tv, 
# 	v_min, v_max, 
# 	v_vals, 
# 	v_tick_labels, 
# 	plot_times, X, 
# 	title, file_name
# 	)

# 	axis = create_axis(title, v_min, v_max, v_vals, v_tick_labels, raw"v")
	
# 	#styles = ["loosely dotted", "loosely dashed", "loosely dashdotted", "loosely dashdotdotted", "solid"]
# 	colors = ["red", "red!80!blue", "red!60!blue", "red!40!blue", "red!20!blue", "blue"]
# 	marks = ["|", "Mercedes star", "x", "star", "asterisk"]
# 	legends = [raw"t = "*string(t) for t = plot_times]
# 	i = 1
	
# 	for t = plot_times
# 		p = tv(t + 0.0025)[1,:]
# 		push!(axis, @pgf PGFPlotsX.Plot(
# 			{
# 				color = colors[i],
# 				#"only marks",
# 				mark = marks[i],
# 				mark_repeat = 6,
# 				mark_options = raw"{scale = 0.9}",
# 				x = "x",
# 				y = "v",
# 			},
# 			Table([:x => X, :v => p])
# 			)
# 		)
# 		i = i + 1
# 	end

# 	push!(axis, @pgf PGFPlotsX.Legend(legends))
	
# 	pgfsave("graphs/pdf/"*file_name*".pdf", axis)
# 	pgfsave("graphs/tikz/"*file_name*".tikz", axis)
	
# end

# function pgfplot_convergence_error(
# 	errors,
# 	mesh_sizes,
# 	error_vals,
# 	error_labels,
# 	title
# )
# 	log_er = log10.(errors.^(-1))
# 	log_ms = log10.(mesh_sizes.^(-1))

# 	(log_er_max, log_er_min) = extrema(log_er)
# 	(log_ms_max, log_ms_min) = extrema(log_ms)

# 	axis = create_axis(
# 		title, 
# 		log_er_max, log_er_min,
# 		error_vals, error_labels,
# 		raw"log_{10}(|u_{h} - u_{{ex}}|_{L^p})"; 
# 		x_min = log_ms_min, x_max = log_ms_min)
	
# 	push!(axis, @pgf PGFPlotsX.Plot(
# 		{
# 			color = "black",
# 			#"only marks",
# 			mark = "x",
# 			#mark_options = raw"{scale = 0.9}",
# 			x = "h",
# 			y = "e",
# 		},
# 		Table([:h => log_ms, :e => log_er])
# 		)
# 	)
# 	pgfsave("graphs/pdf/"*file_name*".pdf", axis)
# 	pgfsave("graphs/tikz/"*file_name*".tikz", axis)
	
# end

# function create_xt(
# 	title, 
# 	t_min, t_max, 
# 	t_vals, 
# 	t_tick_labels
# 	)
# 	"""
# 		Create stereotypical x-t plots for the lambda variable
# 	"""
# 	ytick_str = replace(string(t_vals), "[" => "", "]" => "")

# 	@pgf axis = PGFPlotsX.Axis(
# 		{
# 		view = (0, 90),
# 		colorbar, 
# 		"colormap/jet",
# 		height = raw"12cm",
# 		width = raw"7cm",
# 		title = title,
# 		legend_style = raw"{legend pos = outer north east}",
# 		tick_label_style=raw"{font=small}",
	
# 		xmin = 0,
# 		xmax = 1,
# 		xtick = raw"{0,0.5,1}",
# 		xticklabels = raw"{0,x_0 = frac{1}{2}, 1}",
# 		xlabel_style = raw"{at = {(xticklabel cs:1,0)}}",
	
# 		ymin = t_min,
# 		ymax = t_max,
# 		ytick = raw"{"*ytick_str*"}",
# 		yticklabels = t_tick_labels,
	
# 		xlabel = raw"x",
# 		ylabel = raw"t",
# 		ylabel_style = raw"{at = {(yticklabel cs:1,0)}, rotate = -90}",
	
# 		axis_x_line = "bottom",
# 		axis_y_line = "left",
# 		}
# 	)

# 	axis
# end

# function pgfplot_lambda_1D(
# 	tlambda, 
# 	log10_plot_times, 
# 	X, h, 
# 	title, file_name
# 	)

# 	t_min = minimum(log10_plot_times)
# 	t_max = maximum(log10_plot_times)
# 	axis = create_xt(title, t_min, t_max, [-3, -2, -1, 0], raw"{10^{-3}, 10^{-2}, 10^{-1}, 1}")

# 	function map_tlambda(t,x)
# 		tlambda(t)[Int(round(x/h)) + 1]
# 	end

# 	plot_times = (10).^(log10_plot_times)

# 	push!(axis, @pgf Plot3(
# 		{
# 		surf, 
# 		shader = "flat", 
# 		}, 
# 	PGFPlotsX.Coordinates(X, log10_plot_times, map_tlambda.(plot_times', X))
# 	))

# 	pgfsave("graphs/pdf/"*file_name*".pdf", axis)
# 	pgfsave("graphs/tikz/"*file_name*".tikz", axis)
	
# end