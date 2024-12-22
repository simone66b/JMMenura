using Pkg
Pkg.activate(".")

# Set to location of JMMenura
include(#"./../../../JMMenura/src/JMMenura.jl")

# Import packages
using .JMMenura
using Phylo, Distributions, Random, JLD2, LinearAlgebra, PosDefManifold, Distances, StatsBase, GpABC

#####################
# Set up parameters #
#####################

# Number of traits to simulate
n = # SET

# Open tree - Set file path to tree
tree1 = open(parsenewick, #"./../..//anoles_data//bigsim.tre")

# Find the root node of the trait
root = getroot(tree1)
root_num = tree1.nodedict[root.name]

# Define G matrix - Can be done beforehand and then loaded in using JLD2
@load #"./../../anoles_data/P0.jld2"

# Set trait parameters needed to evolve traits
# Should be a vector of length n 
trait_alpha = #repeat([0.0], n)
trait_mu = #repeat([0.0], n)
trait_sigma = #repeat([sqrt(2)], n)

# Create trait dictionary
# Can allow for different traits to be used on different branches
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

############################
# Set up OU reference data #
############################

# Variables needed for OU matrix model
mat_alpha = #0
mat_sigma = #sqrt(2)
mat_mu = #copy(P0)

# Create matrix dictionary
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

# Specify what evolution functions are used
# Can set dt the evolution lenght
trait_evol_func = trait_evol(#dt = 0.005)
mat_evol_func = mat_evol_affine(#dt = 0.005)

# Run reference simulation
ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func
, mat_evol_func, 0.0, trait_mu, P0, true)

# Extract and save data
ref_data = get_data(ref_sim)
@save "ref_data.jld2" ref_data
@save "tree.jld2" tree1

#####################
# Set up parameters #
#####################

prior = Gamma(2, 0.25)
para = JMMABCAlphaDifferentConstant([prior for _ in 1:n], trait_mu, trait_sigma, prior, mat_mu, mat_sigma, n)

######################
# Perform simulation #
######################

# Set number of points to generate threshold for
thresholds = test_threshold(para_ref_data, tree1, para, zeros(n), P0, 1000, dt = 0.005, each = true)
threshold = sort(thresholds)[4]

# Save threshold
@save "./threshold.jld2" threshold
s_para_result = Vector{Any}()
@save "./s_para_result.jld2" s_para_result