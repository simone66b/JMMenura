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

# Define G matrix - Can be done before hand and then loaded in using JLD2
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

# Load reference data
@load "ref_data.jld2"

#################
# Create Models #
#################
OU_prior = Gamma(2, 0.25)
OU_para = JMMABCAlphaEqualConstant(OU_prior, trait_mu, trait_sigma, OU_prior, mat_mu, mat_sigma, n)

BW_prior = Gamma(2, 0.25)
BW_para = JMMABCBrownian(BW_prior, trait_mu, trait_sigma, mat_mu, BW_prior, n)

Iso_prior = Gamma(2, 0.25)
Iso_para = JMMABCIsospectralAlphaAB(Iso_prior, trait_mu, trait_sigma, Iso_prior, Iso_prior, n)

# Creates functions used to simulate different models
OU_func = create_bayesian_sim(tree1, OU_para, trait_mu, mat_mu, dt = 0.005, each = true)
BW_func = create_bayesian_sim(tree1, BW_para, trait_mu, mat_mu, dt = 0.005, each = true)
Iso_func = create_bayesian_sim(tree1, Iso_para, trait_mu, mat_mu, dt = 0.005, each = true)

model_sim_functions = [OU_func, BW_func, Iso_func]
priors = [get_priors(OU_para), get_priors(BW_para), get_priors(Iso_para)]

# Define function which defines distance inputs.
dist_func = trait_mat_distance(n, nleaves(tree1))

###########################
# Perform Model Selection #
###########################

# Set number of particles to accept and load threshold
n_particles = 200
@load "threshold.jld2"

# Simulate number of points
sim_out = SimulatedModelSelection(ref_data, model_sim_functions, priors, [threshold], n_particles, distance_function = dist_func,
max_iter = 20000*n_particles)

# Save model results
@load "./s_model_result.jld2"
push!(s_model_result, sim_out)
@save "s_model_result.jld2" s_model_result
