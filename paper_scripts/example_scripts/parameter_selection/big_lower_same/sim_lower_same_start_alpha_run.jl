using Pkg

Pkg.activate(".")

include("./../../../../../JMMenura/src/JMMenura.jl")
using .JMMenura
using Phylo, Distributions, Random, JLD2, LinearAlgebra, PosDefManifold, Distances, StatsBase

# number of parameters
n = 4

# Creating tree
tree1 = open(parsenewick, "./../../../..//anoles_data//bigsim.tre")

time_tot = 1.0
tspan = (0.0, time_tot)

# Get root number
root = getroot(tree1)
root_num = tree1.nodedict[root.name]

# G matrix
@load "./../../../../anoles_data/P0.jld2"

@load "./../../../../anoles_data/P1.jld2"

# traits needed to evolve traits
alpha1 = [0.25, 0.5, 1.0, 2.0]
mu1 = repeat([0.0], n)
sigma1 = repeat([sqrt(2)], n)

# create trait dictionary
trait_parameters_true = Dict(root_num => (alpha = alpha1, mu = mu1, sigma = sigma1))
trait_parameters = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha = 1
mat_sigma = sqrt(2)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

@load "para_ref_data.jld2"

#####################
# Set up parameters #
#####################

prior = Gamma(2, 0.25)
para = JMMABCAlphaDifferentConstant([prior for _ in 1:n], mu1, sigma1, prior, mat_mu, mat_sigma, n)

###############
# Warm up sim #
###############

@load "./threshold.jld2" threshold

n_particles = 100

run_result = menura_bayesian(para_ref_data, tree1, para, zeros(n), P0, threshold, n_particles, dt = 0.005, max_iter = 50000*n_particles, each = true)

@load "./s_big_lower_same_result.jld2" sim_lower_same_result

push!(sim_lower_same_result, run_result)

@save "s_big_lower_same_result.jld2" sim_lower_same_result

