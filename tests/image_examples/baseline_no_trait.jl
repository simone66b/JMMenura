# Creates a basic evolution plot

include("../../src/JMMenura.jl")
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold

pyplot()

Random.seed!(1)

# number of parameters
n = 2

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(1, Matrix(1I, n, n)  )))

# traits needed to evolve traits
alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n)
sigma1 = repeat([1.0], n)

# create trait dictionary
trait_parameters = Dict(11 => (alpha = alpha1, mu = mu1, sigma = sigma1))

# Variables needed for OU matrix model
mat_alpha = 1 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

trait_evol_func, mat_evol_func = no_trait_evol()

menura_parameter_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)