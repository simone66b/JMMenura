## cd("/Users/coope/OneDrive/Documents/Uni/Phylogenetics_coding/importance_scripts")
cd("/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/ISO_larger_trait_3")

## include("./../../../JMMenura/src/JMMenura.jl")
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions
using JMMenura


function get_data2(sim_data)
    tree = sim_data[1]
    traits = [tip.data["trait_trace"][end] for tip in getleaves(tree)]
    mats = [tip.data["mat_trace"][end] for tip in getleaves(tree)]
    data = [traits..., mats...]
    return data
end

##########################
# Create Simulation Data #
##########################

n = 4

# Creating tree
trees = open(parsenexus, "fiveTrees.tre")
## tree1 = open(parsenewick, "./anoles_data/bigsim.tre")
time_tot = 1.0
tspan = (0.0, time_tot)
for i in keys(trees) 
# Get root number
root = getroot(trees[i])
root_num = trees[i].nodedict[root.name]
# G matrix
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P1.jld2"

# traits needed to evolve traits
trait_alpha = repeat([0.0], n)
trait_mu = repeat([0.0], n)
trait_sigma = [2, 4, 6, 8]

# create trait dictionary
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))
trait_parameters = (mu = trait_mu, sigma = trait_sigma)

# Variables needed for OU matrix model
mat_a = 5
mat_b = 5

# create matrix dictionary
mat_parameters_true = Dict(root_num => (a = mat_a, b = mat_b))
##mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol_isospectral(dt = 0.01)
trait_evol_func = trait_evol(dt = 0.01)

# Same starting conditions as OU
start_trait_alpha = [2, 4, 6, 8]
start_trait_mu = repeat([0.0], n)
start_trait_sigma = repeat([sqrt(2)], n)

trait_start = start_trait_mu + 3*(start_trait_sigma./sqrt.(2*start_trait_alpha))
mat_start = P1

ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, trees[i], trait_evol_func, mat_evol_func, 0.0, trait_start, P1, true)

## para_ref_data = get_data2(ther_ref_sim)
para_ref_data_tree = push!(para_ref_data_tree, [trees[i], get_data2(ther_ref_sim)])
end ## for loop
@save "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/ISO_larger_trait_3/ISO_larger_trait_3_para_ref_data5D.jld2" para_ref_data_tree
