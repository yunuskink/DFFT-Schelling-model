using StatsBase, StatsPlots, SpecialFunctions, Plots, Random

## OVERVIEW
# Three sims with same segregation index Binary (No DFFT yet)
# Three sims with same segregation index Trinary with vacancies (No DFFT yet)
# Three sims with same segregation index Trinary with three agents (No DFFT yet)
# Intro to Headache function Binary
# Intro to Headache function Trinary with vacancies
# Intro to Headache function Trinary with three agents
# Interpretation of Headache function binary
# Interpretation of Headache function trinary with vacancies
# Interpretation of Headache function trinary with three agents
# Compositional invariance of Headache function binary
# Compositional invariance of Headache function trinary with vacancies
# Compositional invariance of Headache function trinary with three agents
# Forecasting dynamics with binary (transition matrix)
# Forecasting dynamics with trinary (evolution of the mean)

## Intro to Schelling models

##############################
#USER DEFINED PARAMETERS
###########BINARY##############
frac_blue_binary = 0.5
utility_function_blue_binary = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_binary = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

###########BINARY WITH vacancies##############
frac_blue_binary_w_vacancies = 0.33
frac_red_binary_w_vacancies = 0.33
utility_function_blue_binary_w_vacancies = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_binary_w_vacancies = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

###########TRINARY##############
frac_blue_trinary = 0.33
frac_red_trinary = 0.33
utility_function_blue_trinary = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_trinary = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
utility_function_green_trinary = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+2

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_single_binary, lattice_snapshots_binary, utility_function_blue_binary, utility_function_red_binary =
    run_schelling_sim_binary(;
    lattice_length = 120,
    frac_blue_agents = frac_blue_binary,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_binary,
    utility_function_red = utility_function_red_binary,
    bin_length = 5
    )

counts_joint_binary_w_vacancies, lattice_snapshots_binary_w_vacancies,utility_function_blue_binary_w_vacancies,utility_function_red_binary_w_vacancies =
    run_schelling_sim_binary_vacancies(;
    lattice_length = 120,
    frac_red_agents = frac_red_binary_w_vacancies,
    frac_blue_agents = frac_blue_binary_w_vacancies,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_binary_w_vacancies,
    utility_function_red = utility_function_red_binary_w_vacancies,
    bin_length = 5
    )

counts_joint_trinary, lattice_snapshots_trinary,utility_function_blue_trinary,utility_function_red_trinary, utility_function_green_trinary =
    run_schelling_sim_trinary(;
    lattice_length = 120,
    frac_red_agents = frac_red_trinary,
    frac_blue_agents = frac_blue_trinary,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_trinary,
    utility_function_red = utility_function_red_trinary,
    utility_function_green = utility_function_green_trinary,
    bin_length = 5
    )

num_snapshots = size(lattice_snapshots_binary,3)
anim = @animate for i ∈ 1:num_snapshots
    plot_all_three_schelling_sims(
        lattice_snapshots_binary[:,:,i],counts_single_binary,utility_function_blue,utility_function_red,
        lattice_snapshots_binary_w_vacancies[:,:,i],counts_joint_binary_w_vacancies,utility_function_blue_binary_w_vacancies,utility_function_red_binary_w_vacancies,
        lattice_snapshots_trinary[:,:,i],counts_joint_trinary,utility_function_blue_trinary,utility_function_red_trinary,utility_function_green_trinary)
end

gif(anim, "tmp.gif", fps = 2)

## Three sims with same segregation index Binary (No DFFT yet)
# Show three simulations with three different utiity functions that give three different probability distributions but the same segregation indices

##############################
#USER DEFINED PARAMETERS

###########CITY A ##############
utility_function_blue_A = 0.35.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_A = 0.35.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

###########CITY B ##############
utility_function_blue_B = 0.74.*[0, 0, 0, 0, 0, 0, 2, 4, 4]
utility_function_red_B = 0.74.*[0, 0, 0, 0, 0, 0, 2, 4, 4].+ 1

###########CITY C ##############
utility_function_blue_C = 1.0.*[0, 0, 0, 0, 4, 4, 4, 4, 4]
utility_function_red_C = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+ 1

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_single_A, lattice_snapshots_A, utility_function_blue_A, utility_function_red_A =
    run_schelling_sim_binary(;
    lattice_length = 120,
    frac_blue_agents = 0.5,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_A,
    utility_function_red = utility_function_red_A,
    bin_length = 5
    )

counts_single_B, lattice_snapshots_B, utility_function_blue_B, utility_function_red_B =
    run_schelling_sim_binary(;
    lattice_length = 120,
    frac_blue_agents = 0.5,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_B,
    utility_function_red = utility_function_red_B,
    bin_length = 5
    )

counts_single_C, lattice_snapshots_C, utility_function_blue_C, utility_function_red_C =
    run_schelling_sim_binary(;
    lattice_length = 120,
    frac_blue_agents = 0.5,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_C,
    utility_function_red = utility_function_red_C,
    bin_length = 5
    )

num_snapshots = size(lattice_snapshots_A,3)
anim = @animate for i ∈ 1:num_snapshots
    segregation_indices_distributions_comparison_binary_figure(lattice_snapshots_A[:,:,i],counts_single_A,utility_function_blue_A,utility_function_red_A,
        lattice_snapshots_B[:,:,i],counts_single_B,utility_function_blue_B,utility_function_red_B,
        lattice_snapshots_C[:,:,i],counts_single_C,utility_function_blue_C,utility_function_red_C
        )
    # circleplot(x, y, i, line_z = 1:n, cbar = false, c = :reds, framestyle = :none)
end

gif(anim, "tmp.gif", fps = 2)

# segregation_indices_distributions_comparison_binary_figure(lattice_snapshots_A,counts_single_A,utility_function_blue_A,utility_function_red_A,
#     lattice_snapshots_B,counts_single_B,utility_function_blue_B,utility_function_red_B,
#     lattice_snapshots_C,counts_single_C,utility_function_blue_C,utility_function_red_C
#     )

## Three sims with same segregation index Binary with VACANCIES (No DFFT yet)
# Show three simulations with three different utiity functions that give three different probability distributions but the same segregation indices

##############################
#USER DEFINED PARAMETERS

###########CITY A ##############
utility_function_blue_A = 0.5.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_A = 0.5.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

###########CITY B ##############
utility_function_blue_B = 1.5 .*[0, 0, 0, 0.1, 0.25, 1, 4, 9, 9]
utility_function_red_B = 1.5 .*[0, 0, 0, 0.1, 0.25, 1, 4, 9, 9].+ 1

###########CITY C ##############
utility_function_blue_C = 1.0.*[0, 0, 0, 0, 4, 4, 4, 4, 4]
utility_function_red_C = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+ 1

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint_A, lattice_snapshots_A,utility_function_blue_A,utility_function_red_A =
    run_schelling_sim_binary_vacancies(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_A,
    utility_function_red = utility_function_red_A,
    bin_length = 5
    )

counts_joint_B, lattice_snapshots_B,utility_function_blue_B,utility_function_red_B =
    run_schelling_sim_binary_vacancies(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_B,
    utility_function_red = utility_function_red_B,
    bin_length = 5
    )

counts_joint_C, lattice_snapshots_C ,utility_function_blue_C,utility_function_red_C =
    run_schelling_sim_binary_vacancies(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_C,
    utility_function_red = utility_function_red_C,
    bin_length = 5
    )

num_snapshots = size(lattice_snapshots_A,3)
anim = @animate for i ∈ 1:num_snapshots
    segregation_indices_distributions_comparison_binary_w_vacancies(lattice_snapshots_A[:,:,i],counts_joint_A,utility_function_blue_A,utility_function_red_A,
        lattice_snapshots_B[:,:,i],counts_joint_B,utility_function_blue_B,utility_function_red_B,
        lattice_snapshots_C[:,:,i],counts_joint_C,utility_function_blue_C,utility_function_red_C
        )
    # circleplot(x, y, i, line_z = 1:n, cbar = false, c = :reds, framestyle = :none)
end

gif(anim, "tmp.gif", fps = 2)

## Three sims with same segregation index TRINARY (No DFFT yet)
# Show three simulations with three different utiity functions that give three different probability distributions but the same segregation indices

##############################
#USER DEFINED PARAMETERS

###########CITY A ##############
utility_function_blue_A = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red_A = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
utility_function_green_A = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+2

###########CITY B ##############
utility_function_blue_B = 1.5 .*[0, 0, 0, 0.1, 0.25, 1, 4, 9, 9]
utility_function_red_B = 1.5 .*[0, 0, 0, 0.1, 0.25, 1, 4, 9, 9].+ 1
utility_function_green_B = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

###########CITY C ##############
utility_function_blue_C = 1.0.*[0, 0, 0, 0, 4, 4, 4, 4, 4]
utility_function_red_C = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+ 1
utility_function_green_C = 0.5.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint_A, lattice_snapshots_A,utility_function_blue_A,utility_function_red_A =
    run_schelling_sim_trinary(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_A,
    utility_function_red = utility_function_red_A,
    utility_function_green = utility_function_green_A,
    bin_length = 5
    )

counts_joint_B, lattice_snapshots_B,utility_function_blue_B,utility_function_red_B =
    run_schelling_sim_trinary(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_B,
    utility_function_red = utility_function_red_B,
    utility_function_green = utility_function_green_B,
    bin_length = 5
    )

counts_joint_C, lattice_snapshots_C ,utility_function_blue_C,utility_function_red_C =
    run_schelling_sim_trinary(;
    lattice_length = 120,
    frac_red_agents = 0.33,
    frac_blue_agents = 0.33,
    num_simulation_steps = 100,
    utility_function_blue = utility_function_blue_C,
    utility_function_red = utility_function_red_C,
    utility_function_green = utility_function_green_C,
    bin_length = 5
    )

num_snapshots = size(lattice_snapshots_A,3)
anim = @animate for i ∈ 1:num_snapshots
    segregation_indices_distributions_comparison_trinary(
        lattice_snapshots_A[:,:,i],counts_joint_A,utility_function_blue_A,utility_function_red_A,utility_function_green_A,
        lattice_snapshots_B[:,:,i],counts_joint_B,utility_function_blue_B,utility_function_red_B,utility_function_green_B,
        lattice_snapshots_C[:,:,i],counts_joint_C,utility_function_blue_C,utility_function_red_C,utility_function_green_C,
        )
    # circleplot(x, y, i, line_z = 1:n, cbar = false, c = :reds, framestyle = :none)
end
gif(anim, "tmp.gif", fps = 2)

## Intro to Headache function Binary

##############################
#USER DEFINED PARAMETERS
frac_blue = 0.3
utility_function_blue = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_single, lattice_snapshots, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    bin_length = bin_length
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    intro_headache_binary_figures(lattice_snapshots[:,:,i],counts_single,utility_function_blue,utility_function_red)
end
gif(anim, "tmp.gif", fps = 2)

## Intro to Headache function Binary with vacancies
#######################################
#USER DEFINED PARAMETERS
frac_blue = 0.4
frac_red = 0.4
utility_function_blue = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint, lattice_snapshots,utility_function_blue,utility_function_red =
    run_schelling_sim_binary_vacancies(;
    frac_red_agents = frac_red,
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    intro_headache_binary_w_vacancies_figures(lattice_snapshots[:,:,i],counts_joint,utility_function_blue,utility_function_red)
end
gif(anim, "tmp.gif", fps = 2)

## Intro to Headache function Trinary
#######################################
#USER DEFINED PARAMETERS
frac_blue = 0.33
frac_red = 0.33
utility_function_blue = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
utility_function_green = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+2
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint, lattice_snapshots,utility_function_blue,utility_function_red =
    run_schelling_sim_trinary(;
    frac_red_agents = frac_red,
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    utility_function_green = utility_function_green,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    intro_headache_trinary_figures(lattice_snapshots[:,:,i],counts_joint,utility_function_blue,utility_function_red,utility_function_green)
end
gif(anim, "tmp.gif", fps = 2)


## Interpretation of H function Binary

##############################
#USER DEFINED PARAMETERS
#Caution: For now, the forecast in the lower-right panel only works if every possible composition has been observed. So use a small bin_length and segregated utility functions and a large value for t_interval to guarantee this.
frac_blue = 0.4
utility_function_blue = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.3.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length
t_interval = 100
num_simulation_steps = 1000

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_single, lattice_snapshots, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    bin_length = bin_length,
    num_simulation_steps = num_simulation_steps
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    interpretation_headache_binary_figures(lattice_snapshots[:,:,i],counts_single,utility_function_blue,utility_function_red,Int(1+i*t_interval))
end
gif(anim, "tmp.gif", fps = 2)

## Interpretation of H function Binary with Vacancies
#######################################
#USER DEFINED PARAMETERS
frac_blue = 0.33
frac_red = 0.33
utility_function_blue = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.4.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length
t_interval = 2
num_simulation_steps = 100000
initial_number_blue = 2
initial_number_red = 20

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint, lattice_snapshots,utility_function_blue,utility_function_red =
    run_schelling_sim_binary_vacancies(;
    frac_red_agents = frac_red,
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    interpretation_headache_binary_w_vacancies_figures(lattice_snapshots[:,:,i],
        counts_joint,
        utility_function_blue,
        utility_function_red,
        Int(1+i*t_interval),
        initial_number_blue,
        initial_number_red,
        )
end
gif(anim, "tmp.gif", fps = 2)

## Interpretation of H function Trinary
#######################################
#USER DEFINED PARAMETERS
frac_blue = 0.33
frac_red = 0.33
utility_function_blue = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8]
utility_function_red = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+1
utility_function_green = 0.0.*[0, 1, 2, 3, 4, 5, 6, 7, 8].+2
bin_length = 5 #lattice_length of default 60 should be divisible by bin_length
t_interval = 5
num_simulation_steps = 10000
initial_number_blue = 5
initial_number_red = 5

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")

counts_joint, lattice_snapshots,utility_function_blue,utility_function_red =
    run_schelling_sim_trinary(;
    frac_red_agents = frac_red,
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    utility_function_green = utility_function_green,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots,3)
    interpretation_headache_trinary_figures(lattice_snapshots[:,:,i],
        counts_joint,
        utility_function_blue,
        utility_function_red,
        utility_function_green,
        Int(1+i*t_interval),
        initial_number_blue,
        initial_number_red,
        )
end
gif(anim, "tmp.gif", fps = 2)


## Compositional invariance
# An essential trait for any metric of segregation is that it should not depend on...
#1) The overall composition of the region of interest. Otherwise, it would be impossible
#to compare, say, Los Angeles (??% White) with Tucson (??% White). Indeed, any
#information would be confounded by this correlation.
## BINARY Composition Invariance
#######################################
#USER DEFINED PARAMETERS
include("schelling_functions.jl")

fraction_blue_1 = 0.3
fraction_blue_2 = 0.5
fraction_blue_3 = 0.7

utility_function_blue = 0.6 .* [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1
utility_function_red = 0.0 .* [0, 1, 2, 3, 4, 5, 6, 7, 8]

##############################
#SIMULATION AND PLOTTING CODE
counts_single_1, lattice_snapshots_1, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = fraction_blue_1,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

counts_single_2, lattice_snapshots_2, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = fraction_blue_2,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

counts_single_3, lattice_snapshots_3, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = fraction_blue_3,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots_1,3)
    compositional_invariance_binary_figure(lattice_snapshots_1[:,:,i],counts_single_1,fraction_blue_1,
        lattice_snapshots_2[:,:,i],counts_single_2,fraction_blue_2,
        lattice_snapshots_3[:,:,i],counts_single_3,fraction_blue_3,
        utility_function_blue,utility_function_red)
end
gif(anim, "tmp.gif", fps = 2)

## BINARY with VACANCIES Composition Invariance
#######################################
#USER DEFINED PARAMETERS
include("schelling_functions.jl")

utility_function_blue = 0.6 .* [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1
utility_function_red = 0.0 .* [0, 1, 2, 3, 4, 5, 6, 7, 8]

frac_blue_1 = 0.3
frac_red_1 = 0.3

frac_blue_2 = 0.5
frac_red_2 = 0.5

frac_blue_3 = 0.7
frac_red_3 = 0.7

##############################
#SIMULATION AND PLOTTING CODE
counts_joint_1, lattice_snapshots_1,utility_function_blue,utility_function_red =
    run_schelling_sim_binary_vacancies(;
    frac_red_agents = frac_red_1,
    frac_blue_agents = frac_blue_1,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

counts_joint_2, lattice_snapshots_2,utility_function_blue,utility_function_red =
    run_schelling_sim_binary_vacancies(;
    frac_red_agents = frac_red_2,
    frac_blue_agents = frac_blue_2,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

counts_joint_3, lattice_snapshots_3,utility_function_blue,utility_function_red =
    run_schelling_sim_binary_vacancies(;
    frac_red_agents = frac_red_3,
    frac_blue_agents = frac_blue_3,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    )

anim = @animate for i ∈ 1:size(lattice_snapshots_1,3)
    compositional_invariance_binary_figure(lattice_snapshots_1[:,:,i],counts_single_1,fraction_blue_1,
        lattice_snapshots_2[:,:,i],counts_single_2,fraction_blue_2,
        lattice_snapshots_3[:,:,i],counts_single_3,fraction_blue_3,
        utility_function_blue,utility_function_red)
end
gif(anim, "tmp.gif", fps = 2)


## BINARY Sample Size Invariance
#######################################
#USER DEFINED PARAMETERS
lattice_length = 60
bin_length_1 = 4
bin_length_2 = 6
bin_length_3 = 10

frac_blue = 0.5
utility_function_blue = 0.6 .* [0, 1, 2, 3, 4, 5, 6, 7, 8] .+ 1
utility_function_red = 0.0 .* [0, 1, 2, 3, 4, 5, 6, 7, 8]

##############################
#SIMULATION AND PLOTTING CODE
include("schelling_functions.jl")
counts_single_1, lattice_snapshots_1, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    bin_length = bin_length_1,
    lattice_length = lattice_length
    )

counts_single_2, lattice_snapshots_2, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    bin_length = bin_length_2,
    lattice_length = lattice_length
    )

counts_single_3, lattice_snapshots_3, utility_function_blue, utility_function_red =
    run_schelling_sim_binary(;
    frac_blue_agents = frac_blue,
    utility_function_blue = utility_function_blue,
    utility_function_red = utility_function_red,
    bin_length = bin_length_3,
    lattice_length = lattice_length
    )

anim = @animate for i ∈ 1:size(lattice_snapshots_1,3)
    sample_size_invariance_binary_figure(lattice_snapshots_1[:,:,i],counts_single_1,fraction_blue_1,
        lattice_snapshots_2[:,:,i],counts_single_2,fraction_blue_2,
        lattice_snapshots_3[:,:,i],counts_single_3,fraction_blue_3,
        utility_function_blue,utility_function_red)
end
gif(anim, "tmp.gif", fps = 2)
