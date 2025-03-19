using Pkg
Pkg.activate("/home/simoneb/Desktop/JMMenura")

# Set to location of JMMenura
cd("/home/simoneb/Desktop")
include("../../../src/JMMenura.jl")

# Import packages
using .JMMenura
using Phylo, Distributions, Random, JLD2, LinearAlgebra, PosDefManifold, Distances, StatsBase, GpABC

#####################
# Set up parameters #
#####################

# Number of traits to simulate
n = 4 # SET

# Open tree - Set file path to tree
tree1 = open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim.tre")

# Find the root node of the trait
root = getroot(tree1)
root_num = tree1.nodedict[root.name]

# Define G matrix - Can be done before hand and then loaded in using JLD2
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"

# Set trait parameters needed to evolve traits
# Should be a vector of length n 
trait_alpha = repeat([0.0], n)
trait_mu = repeat([0.0], n)
trait_sigma = repeat([sqrt(2)], n)

# Create trait dictionary
# Can allow for different traits to be used on different branches
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

############################
# Set up OU reference data #
############################

# Variables needed for OU matrix model
mat_alpha = 0
mat_sigma = sqrt(2)
mat_mu = copy(P0)

# Create matrix dictionary
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

@load "ref_data.jld2"

#####################
# Set up parameters #
#####################

prior = Gamma(2, 0.25)
para = JMMABCAlphaDifferentConstant([prior for _ in 1:n], trait_mu, trait_sigma, prior, mat_mu, mat_sigma, n)

###############
# Warm up sim #
###############

# Set number of particles to accept and load threshold
n_particles = 100
@load "./threshold.jld2" threshold

run_result = menura_bayesian(ref_data, tree1, para, zeros(n), P0, threshold, n_particles
, dt = 0.005, max_iter = 50000*n_particles, each = true)

# Save model results
@load "./s_para_result.jld2" s_para_result
push!(s_para_result, run_result)
@save "s_para_result.jld2" s_para_result