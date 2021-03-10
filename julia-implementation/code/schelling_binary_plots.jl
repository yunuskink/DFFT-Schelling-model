using Plots
using Plots.PlotMeasures

function utility_function_plot(utility_function_blue,utility_function_red;utility_function_green = [],legend_choice = :topleft)
    # trace_blue = scatter(;x=0:8, y=utility_function_blue, mode="lines+markers")
    # trace_red = scatter(;x=0:8, y=utility_function_red, mode="lines+markers")
    if isempty(utility_function_green)
        plot([0:8],[utility_function_blue,utility_function_red],
            label = ["Blue" "Red"],
            linewidth = 2,
            linecolor = [:blue :red],
            legend = legend_choice,
        )
    else
        plot([0:8],[utility_function_blue,utility_function_red,utility_function_green],
            label = ["Blue" "Red" "Green"],
            linewidth = 2,
            linecolor = [:blue :red :green],
            legend = legend_choice,
        )
    end
end

function lattice_plot_binary(lattice,counts_single;plot_grid = false)
    lattice_length = size(lattice,1)
    p = plot([0.5;lattice_length+0.5;lattice_length+0.5;0.5;0.5],[0.5;0.5;lattice_length+0.5;lattice_length+0.5;0.5],
        framestyle = :none,
        aspect_ratio = 1,
        line = :black,
        hover = false,
        legend = false
        )
    heatmap!(p,lattice,c = :bluesreds)
    if plot_grid
        bin_length = Int(sqrt(length(counts_single)-1))
        num_lines = Int((lattice_length/bin_length) - 1)
        for i=1:num_lines
            coord = (i)*bin_length+0.5
            plot!(p,[0, lattice_length+1],[coord, coord],
            linecolor = :black,
            linewidth = 0.5,
            hover = false,
            alpha = 1,
            legend = false
            )
            plot!(p,[coord, coord],[0, lattice_length+1],
            linecolor = :black,
            linewidth = 0.5,
            hover = false,
            alpha = 1,
            legend = false
            )
        end
    end
    # blue_locs = findall(x->x==1,lattice)
    # red_locs = findall(x->x==2,lattice)
    #
    # sz = 0.4
    # blue_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in blue_locs]
    # blue_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in blue_locs]
    # blue_vertices_x = reduce(vcat,blue_vertices_x)
    # blue_vertices_y = reduce(vcat,blue_vertices_y)
    # plot!(blue_vertices_x,blue_vertices_y,seriestype = :shape,
    #     hover = nothing,
    #     fillcolor = :blue,
    #     linewidth = 0,
    #     legend = false
    #     )
    #
    # red_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in red_locs]
    # red_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in red_locs]
    # red_vertices_x = reduce(vcat,red_vertices_x)
    # red_vertices_y = reduce(vcat,red_vertices_y)
    # plot!(red_vertices_x,red_vertices_y,seriestype = :shape,
    #     hover = nothing,
    #     fillcolor = :red,
    #     linewidth = 0,
    #     legend = false
    #     )
    p
end

function lattice_plot_binary_archive(lattice,counts_single;plot_grid = false)
    lattice_length = size(lattice,1)
    p = plot([0.5;lattice_length+0.5;lattice_length+0.5;0.5;0.5],[0.5;0.5;lattice_length+0.5;lattice_length+0.5;0.5],
        framestyle = :none,
        aspect_ratio = 1,
        line = :black,
        hover = false,
        legend = false
        )
    if plot_grid
        bin_length = Int(sqrt(length(counts_single)-1))
        num_lines = Int((lattice_length/bin_length) - 1)
        for i=1:num_lines
            coord = (i)*bin_length+0.5
            plot!(p,[0, lattice_length],[coord, coord],
                linecolor = :black,
                linewidth = 2,
                hover = false,
                alpha = 0.7,
                legend = false
                )
            plot!(p,[coord, coord],[0, lattice_length],
                linecolor = :black,
                linewidth = 2,
                hover = false,
                alpha = 0.7,
                legend = false
                )
        end
    end

    blue_locs = findall(x->x==1,lattice)
    red_locs = findall(x->x==2,lattice)

    sz = 0.4
    blue_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in blue_locs]
    blue_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in blue_locs]
    blue_vertices_x = reduce(vcat,blue_vertices_x)
    blue_vertices_y = reduce(vcat,blue_vertices_y)
    plot!(blue_vertices_x,blue_vertices_y,seriestype = :shape,
        hover = nothing,
        fillcolor = :blue,
        linewidth = 0,
        legend = false
        )

    red_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in red_locs]
    red_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in red_locs]
    red_vertices_x = reduce(vcat,red_vertices_x)
    red_vertices_y = reduce(vcat,red_vertices_y)
    plot!(red_vertices_x,red_vertices_y,seriestype = :shape,
        hover = nothing,
        fillcolor = :red,
        linewidth = 0,
        legend = false
        )
    p
end


function lattice_composition_binary(lattice,counts_single;plot_grid = false,
    cmap = nothing, colorbar_choice = false)
    N_max = length(counts_single)-1
    lattice = convert.(Int64,lattice)
    bin_length = Int(sqrt(N_max))
    lattice_length = size(lattice,1)
    bins, bin_IDs = initialize_bins(lattice,lattice_length,bin_length)
    lattice_composition = copy(bins)
    for bin_ID in bin_IDs
        lattice_composition[bins.==bin_ID] .= sum((bins.==bin_ID).&(lattice.==1))
    end
    # BORDER LINE
    p = plot([0.5;lattice_length+0.5;lattice_length+0.5;0.5;0.5],[0.5;0.5;lattice_length+0.5;lattice_length+0.5;0.5],
        framestyle = :none,
        aspect_ratio = 1,
        line = :black,
        hover = false,
        legend = false,
        )
    if isnothing(cmap)
        cmap = :redsblues
    end
    heatmap!(p,lattice_composition, c = cmap,
            colorbar = colorbar_choice,
            colorbar_title = "# Blue"
            )
    # if plot_grid
    #     bin_length = Int(sqrt(length(counts_single)-1))
    #     num_lines = Int((lattice_length/bin_length) - 1)
    #     for i=1:num_lines
    #         coord = (i)*bin_length+0.5
    #         plot!(p,[0, lattice_length],[coord, coord],
    #             linecolor = :black,
    #             linewidth = 2,
    #             hover = false,
    #             alpha = 0.7,
    #             legend = false
    #             )
    #         plot!(p,[coord, coord],[0, lattice_length],
    #             linecolor = :black,
    #             linewidth = 2,
    #             hover = false,
    #             alpha = 0.7,
    #             legend = false
    #             )
    #     end
    # end
    #
    # blue_locs = findall(x->x==1,lattice)
    # red_locs = findall(x->x==2,lattice)
    #
    # sz = 0.4
    # blue_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in blue_locs]
    # blue_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in blue_locs]
    # blue_vertices_x = reduce(vcat,blue_vertices_x)
    # blue_vertices_y = reduce(vcat,blue_vertices_y)
    # plot!(blue_vertices_x,blue_vertices_y,seriestype = :shape,
    #     hover = nothing,
    #     fillcolor = :blue,
    #     linewidth = 0,
    #     legend = false
    #     )
    #
    # red_vertices_x = [[(loc[1]-sz), (loc[1]+sz), (loc[1]+sz),(loc[1]-sz),NaN]  for loc in red_locs]
    # red_vertices_y = [[(loc[2]-sz), (loc[2]-sz), (loc[2]+sz),(loc[2]+sz),NaN]  for loc in red_locs]
    # red_vertices_x = reduce(vcat,red_vertices_x)
    # red_vertices_y = reduce(vcat,red_vertices_y)
    # plot!(red_vertices_x,red_vertices_y,seriestype = :shape,
    #     hover = nothing,
    #     fillcolor = :red,
    #     linewidth = 0,
    #     legend = false
    #     )
    # p
end


function headache_plot_binary(counts_single;
    x_label = "# blue agents",
    linecolor = :blue,
    legend_labels = [""]
    )
    #counts_single is an array with each column representing a different experiment
    H = Array{Float64,2}(undef,length(counts_single),0)
    V_n = Array{Float64,2}(undef,length(counts_single),0)
    f_n = Array{Float64,2}(undef,length(counts_single),0)
    N_max = length(counts_single)-1
    n_s = 0:N_max
    n_s = collect(n_s)./N_max
    H_tmp = calculate_H_binary(counts_single)
    y_int, V_i = fit_linear(n_s,H_tmp)
    H_tmp = H_tmp .- y_int
    y_int, V_i = fit_linear(n_s,H_tmp)
    H = hcat(H,H_tmp)
    V_tmp = V_i.*n_s
    V_n = hcat(V_n,V_tmp)
    f_tmp = H_tmp.-V_tmp
    f_n = hcat(f_n,f_tmp)

    p = plot(H,label = "H",linecolor = linecolor,
        linewidth = 2)
    plot!(p,f_n,linestyle = :dash,label = "f",linecolor = linecolor,
        linewidth = 2)
    plot!(p,V_n,linestyle = :dot,label = "Vn",linecolor = linecolor,
        linewidth = 2)

    xlabel!(p,x_label,legend_choice = :topleft)
    p
end

function rate_plot_binary(counts_single;
    x_label = "# blue agents",
    linecolor = :blue,
    legend_labels = [""]
    )
    #counts_single is an array with each column representing a different experiment
    H = calculate_H_binary(counts_single)
    N_max = length(counts_single)-1
    delta_H = N_max.*(H[1:end-1] .- H[2:end])
    R = 1.0./(1.0.+exp.(delta_H))
    p = plot(R,label = "b->r",
        linecolor = linecolor,
        linewidth = 2)
    xlabel!(p,x_label)
    p
end

function transition_matrix_plot_binary(counts_single,
    t
    )
    #counts_single is an array with each column representing a different experiment
    H = calculate_H_binary(counts_single)
    N_max = length(counts_single)-1
    dH_dN_increasing = N_max.*(H[2:end] .- H[1:end-1])
    dH_dN_decreasing = N_max.*(H[1:end-1] .- H[2:end])
    n_s = collect(0:N_max)./N_max
    decreasing_N_rate = n_s[2:end].*exp.(-dH_dN_decreasing)./(1.0.+exp.(-dH_dN_decreasing));
    increasing_N_rate = (1.0.-n_s[1:end-1]).*exp.(-dH_dN_increasing)./(1.0.+exp.(-dH_dN_increasing));
    constant_N_rate = 1.0 .- [0;decreasing_N_rate] .- [increasing_N_rate;0];
    T_matrix = zeros(size(H,1),size(H,1))
    for i=1:length(decreasing_N_rate)
        T_matrix[i+1,i] = decreasing_N_rate[i]
    end
    for i=1:length(increasing_N_rate)
        T_matrix[i,i+1] = increasing_N_rate[i]
    end
    for i=1:length(constant_N_rate)
        T_matrix[i,i] = constant_N_rate[i]
    end
    # T_matrix = diag(constant_N_rate) + diag(increasing_N_rate,1) + diag(decreasing_N_rate,-1);
    # println(increasing_N_rate)
    # println(constant_N_rate)

    ns_label = [string(i) for i = 0:N_max]
    t = Int(t)
    T_matrix = T_matrix^t
    T_matrix = rotl90(T_matrix)
    T_matrix = reverse(T_matrix,dims = 1)
    # T_matrix = reverse(T_matrix,dims = 1)
    # T_matrix = reverse(T_matrix,dims = 2)
    p = heatmap(ns_label,ns_label,T_matrix,aspect_ratio = 1)

    # p = plot(R,label = "H",
    #     linecolor = linecolor,
    #     linewidth = 2)
    # xlabel!(p,x_label)
    p
end


function segregation_indices_plot_binary(counts_single)
    H_theil, D = calculate_seg_indices_binary(counts_single)
    bar(["Entropy","Dissimilarity"], [H_theil,D],
        legend = false
        )
    # plot(bar_plot)
end

function intro_headache_binary_figures(lattice_snapshots,
    counts_single,
    utility_function_blue,
    utility_function_red)
    # plotly()
    p1 = plot([0:8],[utility_function_blue,utility_function_red],
        label = ["Blue" "Red"],
        linewidth = 2,
        linecolor = [:blue :red],
        legend = :topleft)
    xlabel!(p1,"# neighbors")
    # p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = lattice_composition_binary(lattice_snapshots[:,:,end],counts_single)
    p4 = bar([0:(length(counts_single)-1)],counts_single,legend = false)
    xlabel!(p4,"# blue agents")
    p5 = headache_plot_binary(counts_single)

    plot(p1,p2,p3,p4,p5,
        layout = (1,5),
        title = ["Utility" "Snapshot" "Compositions" "Histogram" "DFFT"],
        size = (900,200),
        bottom_margin = 5mm
    )
end

function interpretation_headache_binary_figures(lattice_snapshots,
    counts_single,
    utility_function_blue,
    utility_function_red,
    t)
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = bar([0:(length(counts_single)-1)],counts_single,legend = false)
    p4 = headache_plot_binary(counts_single)
    p5 = rate_plot_binary(counts_single)
    p6 = transition_matrix_plot_binary(counts_single,t)

    plot(p1,p2,p3,p4,p5,p6,
        layout = (2,3),
        title = ["Utility" "Snapshot" "Histogram" "DFFT" "Rate b->r" string(t)],
    )
end

function dynamics_binary_figures(lattice_snapshots,
    counts_single,
    utility_function_blue,
    utility_function_red,
    T_matrix_schelling,
    T_matrix_DFFT,
    schelling_steps,
    DFFT_steps,
    )
    p1 = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    p2 = lattice_plot_binary(lattice_snapshots[:,:,end],counts_single)
    p3 = bar([0:(length(counts_single)-1)],counts_single,legend = false)
    p4 = headache_plot_binary(counts_single)
    p5 = rate_plot_binary(counts_single)
    N_max = length(counts_single)-1
    # p6 = transition_matrix_plot_binary(counts_single,t)
    ns_label = [string(i) for i = 0:N_max]
    max_prob = maximum([maximum(T_matrix_DFFT);maximum(T_matrix_schelling)])
    p6 = heatmap(ns_label,ns_label,T_matrix_DFFT,aspect_ratio = 1,clim = (0,max_prob))
    p7 = heatmap(ns_label,ns_label,T_matrix_schelling,aspect_ratio = 1,clim = (0,max_prob))
    # l = @layout[a{0.2w} grid(3,2) grid(2,1)]
    l = @layout[grid(1,5); grid(1,2)]
    plot(p1,p2,p3,p4,p5,p6,p7,
        layout = l,
        title = ["Utility" "Snapshot" "Histogram" "DFFT" "Rates" string(DFFT_steps) string(schelling_steps)],
        size = (900,500)
    )
end


function compositional_invariance_binary_figure(lattice_snapshots_1,counts_single_1,fraction_blue_1,
    lattice_snapshots_2,counts_single_2,fraction_blue_2,
    lattice_snapshots_3,counts_single_3,fraction_blue_3,
    utility_function_blue,utility_function_red;
    linecolor_1 = :blue,
    linecolor_2 = :green,
    linecolor_3 = :red,
    )
    # pyplot()

    p_utility = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    xlabel!(p_utility,"# Neighbors")
    # ylabel!(p_utility_A,"A")
    # ylabel!(p_utility_B,"B")
    # ylabel!(p_utility_C,"C")

    p_lattice_1 = lattice_plot_binary(lattice_snapshots_1[:,:,end],counts_single_1;plot_grid = false)
    p_lattice_2 = lattice_plot_binary(lattice_snapshots_2[:,:,end],counts_single_2;plot_grid = false)
    p_lattice_3 = lattice_plot_binary(lattice_snapshots_3[:,:,end],counts_single_3;plot_grid = false)

    p_lattice_comp_1 = lattice_composition_binary(lattice_snapshots_1[:,:,end],counts_single_1;plot_grid = true)
    p_lattice_comp_2 = lattice_composition_binary(lattice_snapshots_2[:,:,end],counts_single_2;plot_grid = true)
    p_lattice_comp_3 = lattice_composition_binary(lattice_snapshots_3[:,:,end],counts_single_3;plot_grid = true)

    # p_distribution_1 = bar([0:(length(counts_single_1)-1)],counts_single_1,legend = false)
    # p_distribution_2 = bar([0:(length(counts_single_2)-1)],counts_single_2,legend = false)
    # p_distribution_3 = bar([0:(length(counts_single_3)-1)],counts_single_3,legend = false)

    p_distribution_1 = bar([0:1/(length(counts_single_1)-1):1],counts_single_1,legend = false,linecolor = :match,xticks = 0:1,label = "Simulation",legendfontsize=7)
    p_distribution_2 = bar([0:1/(length(counts_single_2)-1):1],counts_single_2,legend = false,linecolor = :match,xticks = 0:1,label = "Simulation")
    p_distribution_3 = bar([0:1/(length(counts_single_3)-1):1],counts_single_3,legend = false,linecolor = :match,xticks = 0:1,label = "Simulation")

    # Now plot the naive binary distribution

    N_max_1 = length(counts_single_1)-1
    N_max_2 = length(counts_single_2)-1
    N_max_3 = length(counts_single_3)-1
    Ns_1 = 0:1:N_max_1
    binary_dist_1 = sum(counts_single_1).*[factorial(big(N_max_1))/(factorial(big(N))*factorial(big(N_max_1-N)))*fraction_blue_1^N*(1.0-fraction_blue_1)^(N_max_1-N) for N in Ns_1]
    Ns_2 = 0:1:N_max_2
    binary_dist_2 = sum(counts_single_2).*[factorial(big(N_max_2))/(factorial(big(N))*factorial(big(N_max_2-N)))*fraction_blue_2^N*(1.0-fraction_blue_2)^(N_max_2-N) for N in Ns_2]
    Ns_3 = 0:1:N_max_3
    binary_dist_3 = sum(counts_single_3).*[factorial(big(N_max_3))/(factorial(big(N))*factorial(big(N_max_3-N)))*fraction_blue_3^N*(1.0-fraction_blue_3)^(N_max_3-N) for N in Ns_3]

    plot!(p_distribution_1,[0:1/(length(counts_single_1)-1):1],binary_dist_1,label = "Random")
    plot!(p_distribution_2,[0:1/(length(counts_single_2)-1):1],binary_dist_2,label = "Random")
    plot!(p_distribution_3,[0:1/(length(counts_single_3)-1):1],binary_dist_3,label = "Random")

    # ylim!(p_distribution_1,(0,))

    xlabel!(p_distribution_3,"Fraction Blue")

    H_1 = calculate_H_binary(counts_single_1)
    H_2 = calculate_H_binary(counts_single_2)
    H_3 = calculate_H_binary(counts_single_3)
    # p_headaches_separate = plot([0:(length(H_1)-1)],
    #     [H_1 H_2 H_3],
    #     label = [fraction_blue_1 fraction_blue_2 fraction_blue_3],
    #     )
    # xlabel!(p_headaches_separate,"Fraction Blue")


    p_headaches_separate = plot([0:(length(H_1)-1)]./N_max_1,
        [H_1],
        label = "H(n): "*string(Int(round(100*fraction_blue_1)))*"%",
        linecolor=linecolor_1,
        linewidth = 2,
        linestyle = :dash,
        legend = :outerright,
        legendfontsize = 9
        # legendtitle = "function:composition"
        )
    plot!(p_headaches_separate,
        [0:(length(H_2)-1)]./N_max_2,
        [H_2],
        label = "H(n): "*string(Int(round(100*fraction_blue_2)))*"%",
        linecolor=linecolor_2,
        linewidth = 2,
        linestyle = :dash,
        )
    plot!(p_headaches_separate,
        [0:(length(H_3)-1)]./N_max_3,
        [H_3],
        label = "H(n): "*string(Int(round(100*fraction_blue_3)))*"%",
        linecolor=linecolor_3,
        linewidth = 2,
        linestyle = :dash,
        )


    xlabel!(p_headaches_separate,"Fraction Blue")


    H_arr = [H_1 H_2 H_3]
    H_mean = similar(H_1)
    for i=1:size(H_arr,1)
        H_tmp = 0
        num_not_nan = 0
        for j=1:size(H_arr,2)
            if ~isnan(H_arr[i,j])
                H_tmp += H_arr[i,j]
                num_not_nan += 1
            end
        end
        H_mean[i] = H_tmp./num_not_nan
    end
    N_max = length(H_1)-1
    ns = (0:1/N_max:1)
    C_1,V_1 = fit_linear((0:1/N_max:1),H_1.-H_mean)
    C_2,V_2 = fit_linear((0:1/N_max:1),H_2.-H_mean)
    C_3,V_3 = fit_linear((0:1/N_max:1),H_3.-H_mean)

    f_1 = H_1 .- V_1.*(0:1/N_max:1)
    f_2 = H_2 .- V_2.*(0:1/N_max:1)
    f_3 = H_3 .- V_3.*(0:1/N_max:1)

    C_shift = minimum(f_1[.!isnan.(f_1)])
    f_2_diff = f_2.-f_1
    f_3_diff = f_3.-f_1
    f_1 = f_1 .- C_shift
    f_2 = f_2 .- mean(f_2_diff[.!isnan.(f_2_diff)]).- C_shift
    f_3 = f_3 .- mean(f_3_diff[.!isnan.(f_3_diff)]).- C_shift
    C_1 = C_shift
    C_2 = C_shift + mean(f_2_diff[.!isnan.(f_2_diff)])
    C_3 = C_shift + mean(f_3_diff[.!isnan.(f_3_diff)])


    plot!(p_headaches_separate,
    [0:(length(H_1)-1)]./N_max_1,
    [f_1],
    label = "f(n): "*string(Int(round(100*fraction_blue_1)))*"%",
    linecolor=linecolor_1,
    linewidth = 2,
    linestyle = :dot
    )
    plot!(p_headaches_separate,
    [0:(length(H_2)-1)]./N_max_2,
    [f_2],
    label = "f(n): "*string(Int(round(100*fraction_blue_2)))*"%",
    linecolor=linecolor_2,
    linewidth = 2,
    linestyle = :dot
    )
    plot!(p_headaches_separate,
    [0:(length(H_3)-1)]./N_max_3,
    [f_3],
    label = "f(n): "*string(Int(round(100*fraction_blue_3)))*"%",
    linecolor=linecolor_3,
    linewidth = 2,
    linestyle = :dot
    )

    plot!(p_headaches_separate,
    [0:(length(H_2)-1)]./N_max_2,
    [V_1.*ns.+C_1],
    label =  "Vn+z: "*string(Int(round(100*fraction_blue_1)))*"%",
    linecolor=linecolor_1,
    linewidth = 2,
    linestyle = :solid,
    )
    plot!(p_headaches_separate,
    [0:(length(H_2)-1)]./N_max_2,
    [V_2.*ns.+C_2],
    label = "Vn+z: "*string(Int(round(100*fraction_blue_2)))*"%",
    linecolor=linecolor_2,
    linewidth = 2,
    linestyle = :solid,
    )
    plot!(p_headaches_separate,
    [0:(length(H_3)-1)]./N_max_3,
    [V_3.*ns.+C_3],
    label = "Vn+z: "*string(Int(round(100*fraction_blue_3)))*"%",
    linecolor=linecolor_3,
    linewidth = 2,
    linestyle = :solid,
    )

    # p_DFFT_decomposition = plot(0:N_max,
    #     [f_1 f_2 f_3],
    #     label = [fraction_blue_1 fraction_blue_2 fraction_blue_3]
    #     )
    # xlabel!(p_DFFT_decomposition,"# Blue")
    # plot!(p_DFFT_decomposition,
    #     0:N_max,
    #     [collect(C_1.+V_1.*(0:N_max)) collect(C_2.+V_2.*(0:N_max)) collect(C_3.+V_3.*(0:N_max))],
    #     )
    # p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
    #     bar_position = :dodge,
    #     label = ["Entropy" "Dissimilarity"],
    #     legend = :topright
    #     )
    # ylims!(p_seg_indices_bar,(0,1))
    # xlabel!(p_seg_indices_bar,"City")

    # l = @layout[a{0.2w} grid(3,2) grid(2,1)]
    # plot(p_utility,
    #     p_lattice_1,p_distribution_1,
    #     p_lattice_2,p_distribution_2,
    #     p_lattice_3,p_distribution_3,
    #     p_headaches_separate,
    #     p_DFFT_decomposition,
    #     layout = l,
    #     title = ["Utility" "Snapshot" "Histogram" "" "" "" "" "H(N_b)" "f(N_b)"],
    #     size = (900,500),
    #     bottom_margin = 5mm
    #     )

    l = @layout[a{0.15w} grid(3,3) b{0.3w}]
    plot(p_utility,
        p_lattice_1,p_lattice_comp_1,p_distribution_1,
        p_lattice_2,p_lattice_comp_2,p_distribution_2,
        p_lattice_3,p_lattice_comp_3,p_distribution_3,
        p_headaches_separate,
        layout = l,
        # title = ["Utility" string(bin_length_1,"x",bin_length_1) "Compositions" "Histogram" string(bin_length_2,"x",bin_length_2) "" "" string(bin_length_3,"x",bin_length_3) "" "" "H(n)"],
        title = ["Utility" "Snapshot" "Compositions" "Histogram" "" "" "" "" "" "" "DFFT functions"],
        titlefontsize = 12,
        size = (900,450),
        bottom_margin = 5mm,
        dpi = 600
        )


    #Show the following figures...
    #1) Utility functions
    #2) Snapshots of the three cities
    #3) Three Headache functions
    #4) Three aligned frustration functions and three linear vexation parts
end

function sample_size_invariance_binary_figure(lattice_snapshots_1,counts_single_1,fraction_blue_1,
    lattice_snapshots_2,counts_single_2,fraction_blue_2,
    lattice_snapshots_3,counts_single_3,fraction_blue_3,
    utility_function_blue,utility_function_red)
    # pyplot()

    p_utility = utility_function_plot(utility_function_blue,utility_function_red,legend_choice = false)
    xlabel!(p_utility,"# Neighbors")
    # ylabel!(p_utility_A,"A")
    # ylabel!(p_utility_B,"B")
    # ylabel!(p_utility_C,"C")

    p_lattice_1 = lattice_plot_binary(lattice_snapshots_1[:,:,end],counts_single_1;plot_grid = false)
    p_lattice_2 = lattice_plot_binary(lattice_snapshots_2[:,:,end],counts_single_2;plot_grid = false)
    p_lattice_3 = lattice_plot_binary(lattice_snapshots_3[:,:,end],counts_single_3;plot_grid = false)

    p_lattice_comp_1 = lattice_composition_binary(lattice_snapshots_1[:,:,end],counts_single_1;plot_grid = true)
    p_lattice_comp_2 = lattice_composition_binary(lattice_snapshots_2[:,:,end],counts_single_2;plot_grid = true)
    p_lattice_comp_3 = lattice_composition_binary(lattice_snapshots_3[:,:,end],counts_single_3;plot_grid = true)

    p_distribution_1 = bar([0:1/(length(counts_single_1)-1):1],counts_single_1,legend = false,linecolor = :match,xticks = 0:1)
    p_distribution_2 = bar([0:1/(length(counts_single_2)-1):1],counts_single_2,legend = false,linecolor = :match,xticks = 0:1)
    p_distribution_3 = bar([0:1/(length(counts_single_3)-1):1],counts_single_3,legend = false,linecolor = :match,xticks = 0:1)
    xlabel!(p_distribution_3,"Fraction Blue")

    N_max_1 = length(counts_single_1)-1
    N_max_2 = length(counts_single_2)-1
    N_max_3 = length(counts_single_3)-1
    Ns_1 = 0:1:N_max_1
    binary_dist_1 = sum(counts_single_1).*[factorial(big(N_max_1))/(factorial(big(N))*factorial(big(N_max_1-N)))*fraction_blue_1^N*(1.0-fraction_blue_1)^(N_max_1-N) for N in Ns_1]
    Ns_2 = 0:1:N_max_2
    binary_dist_2 = sum(counts_single_2).*[factorial(big(N_max_2))/(factorial(big(N))*factorial(big(N_max_2-N)))*fraction_blue_2^N*(1.0-fraction_blue_2)^(N_max_2-N) for N in Ns_2]
    Ns_3 = 0:1:N_max_3
    binary_dist_3 = sum(counts_single_3).*[factorial(big(N_max_3))/(factorial(big(N))*factorial(big(N_max_3-N)))*fraction_blue_3^N*(1.0-fraction_blue_3)^(N_max_3-N) for N in Ns_3]

    plot!(p_distribution_1,[0:1/(length(counts_single_1)-1):1],binary_dist_1,label = "Random")
    plot!(p_distribution_2,[0:1/(length(counts_single_2)-1):1],binary_dist_2,label = "Random")
    plot!(p_distribution_3,[0:1/(length(counts_single_3)-1):1],binary_dist_3,label = "Random")


    # p_distribution_1 = bar([0:(length(counts_single_1)-1)],counts_single_1,legend = false)
    # p_distribution_2 = bar([0:(length(counts_single_2)-1)],counts_single_2,legend = false)
    # p_distribution_3 = bar([0:(length(counts_single_3)-1)],counts_single_3,legend = false)
    # xlabel!(p_distribution_3,"# Blue")

    H_1 = calculate_H_binary(counts_single_1)
    H_1 = H_1 .- maximum(H_1[.!isnan.(H_1)])
    H_2 = calculate_H_binary(counts_single_2)
    H_2 = H_2 .- maximum(H_2[.!isnan.(H_2)])
    # H_2_diff = H_1.-H_2
    # H_2 = H_2 .- mean(H_2_diff[.!isnan.(H_2_diff)])
    H_3 = calculate_H_binary(counts_single_3)
    H_3 = H_3 .- maximum(H_3[.!isnan.(H_3)])
    C_shift = minimum([H_1[.!isnan.(H_1)];H_2[.!isnan.(H_2)];H_3[.!isnan.(H_3)]])
    # H_3_diff = H_1.-H_2
    # H_3 = H_3 .- mean(H_3_diff[.!isnan.(H_3_diff)])
    N_max_1 = length(counts_single_1)-1
    N_max_2 = length(counts_single_2)-1
    N_max_3 = length(counts_single_3)-1
    p_headaches_separate = plot([0:(length(H_1)-1)]./N_max_1,
        [H_1.+C_shift],
        label = "H(n)-z: "*string(bin_length_1)*"x"*string(bin_length_1),
        legend = :outertopright,
        legendfontsize = 9,
        linewidth = 3,
        linestyle = :dash,
        )
    plot!(p_headaches_separate,
        [0:(length(H_2)-1)]./N_max_2,
        [H_2.+C_shift],
        label = "H(n)-z: "*string(bin_length_2)*"x"*string(bin_length_2),
        linewidth = 3,
        linestyle = :dash,
        )
    plot!(p_headaches_separate,
        [0:(length(H_3)-1)]./N_max_3,
        [H_3.+C_shift],
        label = "H(n)-z: "*string(bin_length_3)*"x"*string(bin_length_3),
        linewidth = 3,
        linestyle = :dash,
        )

    xlabel!(p_headaches_separate,"Fraction Blue")
    # plot!(p_DFFT_decomposition,
    #     0:N_max,
    #     [collect(C_1.+V_1.*(0:N_max)) collect(C_2.+V_2.*(0:N_max)) collect(C_3.+V_3.*(0:N_max))],
    #     )
    # p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
    #     bar_position = :dodge,
    #     label = ["Entropy" "Dissimilarity"],
    #     legend = :topright
    #     )
    # ylims!(p_seg_indices_bar,(0,1))
    # xlabel!(p_seg_indices_bar,"City")

    l = @layout[a{0.15w} grid(3,3) b{0.3w}]

    plot(p_utility,
        p_lattice_1,p_lattice_comp_1,p_distribution_1,
        p_lattice_2,p_lattice_comp_2,p_distribution_2,
        p_lattice_3,p_lattice_comp_3,p_distribution_3,
        p_headaches_separate,
        layout = l,
        # title = ["Utility" string(bin_length_1,"x",bin_length_1) "Compositions" "Histogram" string(bin_length_2,"x",bin_length_2) "" "" string(bin_length_3,"x",bin_length_3) "" "" "H(n)"],
        title = ["Utility" "Snapshot" "Compositions" "Histogram" "" "" "" "" "" "" "DFFT functions"],
        titlefontsize = 12,
        size = (900,450),
        bottom_margin = 5mm,
        dpi = 600
        )

    #Show the following figures...
    #1) Utility functions
    #2) Snapshots of the three cities
    #3) Three Headache functions
    #4) Three aligned frustration functions and three linear vexation parts
end


function segregation_indices_distributions_comparison_binary_figure(lattice_snapshots_A,counts_single_A,utility_function_blue_A,utility_function_red_A,
    lattice_snapshots_B,counts_single_B,utility_function_blue_B,utility_function_red_B,
    lattice_snapshots_C,counts_single_C,utility_function_blue_C,utility_function_red_C
    )
    # pyplot()
    # gr()
    #Each
    p_utility_A = utility_function_plot(utility_function_blue_A,utility_function_red_A,legend_choice = false)
    p_utility_B = utility_function_plot(utility_function_blue_B,utility_function_red_B,legend_choice = false)
    p_utility_C = utility_function_plot(utility_function_blue_C,utility_function_red_C,legend_choice = false)
    xlabel!(p_utility_C,"# Neighbors")

    p_utility_A = plot(p_utility_A,
        annotation=(-3,(minimum([utility_function_blue_A;utility_function_red_A])+maximum([utility_function_blue_A;utility_function_red_A]))/2,"A"))
    p_utility_B = plot(p_utility_B,
        annotation=(-3,(minimum([utility_function_blue_B;utility_function_red_B])+maximum([utility_function_blue_B;utility_function_red_B]))/2,"B"))
    p_utility_C = plot(p_utility_C,
        annotation=(-3,(minimum([utility_function_blue_C;utility_function_red_C])+maximum([utility_function_blue_C;utility_function_red_C]))/2,"C"))

    # p_utility_A = plot(p_utility_A,
    #     annotation=(-0.1,(minimum([utility_function_blue_A;utility_function_red_A])+maximum([utility_function_blue_A;utility_function_red_A]))/2,"A"))
    # p_utility_A = plot(p_utility_A,
    #     annotation=(-0.1,(minimum([utility_function_blue_A;utility_function_red_A])+maximum([utility_function_blue_A;utility_function_red_A]))/2,"A"))

    # ylabel!(p_utility_A,"A")
    # ylabel!(p_utility_B,"B")
    # ylabel!(p_utility_C,"C")

    p_lattice_A = lattice_plot_binary(lattice_snapshots_A[:,:,end],counts_single_A)
    p_lattice_B = lattice_plot_binary(lattice_snapshots_B[:,:,end],counts_single_B)
    p_lattice_C = lattice_plot_binary(lattice_snapshots_C[:,:,end],counts_single_C)

    p_lattice_comp_A = lattice_composition_binary(lattice_snapshots_A[:,:,end],counts_single_A)
    p_lattice_comp_B = lattice_composition_binary(lattice_snapshots_B[:,:,end],counts_single_B)
    p_lattice_comp_C = lattice_composition_binary(lattice_snapshots_C[:,:,end],counts_single_C)

    p_distribution_A = bar([0:(length(counts_single_C)-1)],counts_single_A,legend = false,linecolor = :match)
    p_distribution_B = bar([0:(length(counts_single_C)-1)],counts_single_B,legend = false,linecolor = :match)
    p_distribution_C = bar([0:(length(counts_single_C)-1)],counts_single_C,legend = false,linecolor = :match)
    xlabel!(p_distribution_C,"# Blue")

    H_theil_A, D_A = calculate_seg_indices_binary(counts_single_A)
    H_theil_B, D_B = calculate_seg_indices_binary(counts_single_B)
    H_theil_C, D_C = calculate_seg_indices_binary(counts_single_C)

    # p_seg_indices_bar = groupedbar(["A";"B";"C"],[H_theil_A D_A;H_theil_B D_B;H_theil_C D_C],
    #     bar_position = :dodge,
    #     label = ["Entropy" "Dissimilarity"],
    #     legend = :topright
    #     )
    # ylims!(p_seg_indices_bar,(0,1))
    # xlabel!(p_seg_indices_bar,"City")

    p_seg_indices_bar_theil = bar(["A";"B";"C"],[H_theil_A;H_theil_B;H_theil_C],
        # bar_position = :dodge,
        # label = ["Entropy" "Dissimilarity"],
        legend = false
        )
    ylims!(p_seg_indices_bar_theil,(0,1))
    ylabel!(p_seg_indices_bar_theil,"Theil segregation index")
    xlabel!(p_seg_indices_bar_theil,"City")

    p_seg_indices_bar_dis = bar(["A";"B";"C"],[D_A;D_B;D_C],
        # bar_position = :dodge,
        # label = ["Entropy" "Dissimilarity"],
        legend = false
        )
    ylims!(p_seg_indices_bar_dis,(0,1))
    ylabel!(p_seg_indices_bar_dis,"Dissimilarity index")
    xlabel!(p_seg_indices_bar_dis,"City")


    l = @layout[grid(3,4) grid(2,1){0.2w}]
    plot(p_utility_A,p_lattice_A,p_lattice_comp_A,p_distribution_A,
        p_utility_B,p_lattice_B,p_lattice_comp_B,p_distribution_B,
        p_utility_C,p_lattice_C,p_lattice_comp_C,p_distribution_C,
        p_seg_indices_bar_dis,p_seg_indices_bar_theil, layout = l,
        title = ["Utility" "Snapshot" "Compositions" "Histogram" "" "" "" "" "" "" "" "" "Seg. index" ""],
        size = (900,400),
        left_margin = 5mm,
        bottom_margin = 5mm
        )
    # l = @layout[grid(3,1){0.15w} grid(3,1){0.2w} grid(3,1){0.3w} grid(3,1){0.15w} grid(2,1){0.2w}]
    # plot(p_utility_A,
    # p_utility_B,
    # p_utility_C,
    # p_lattice_A,
    # p_lattice_B,
    # p_lattice_C,
    # p_lattice_comp_A,
    # p_lattice_comp_B,
    # p_lattice_comp_C,
    # p_distribution_A,
    #     p_distribution_B,
    #     p_distribution_C,
    #     p_seg_indices_bar_theil,
    #     p_seg_indices_bar_dis,
    #      layout = l,
    #     title = ["Utility" "" "" "Snapshot" "" "" "Compositions" "" "" "Histogram" "" "" "Seg. index" ""],
    #     size = (900,400)
    #     )

end
