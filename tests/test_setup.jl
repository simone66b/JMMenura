# Test file which sets up a basic instance for a run
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

# number of parameters
n = 9

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) # Change to Wishart Distribution

# traits needed to evolve traits
alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n) ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n)
trait_para = (alpha=alpha1, sigma=sigma1, mu = mu1)

# Variables needed for OU matrix model
mat_alpha = 0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)
mat_para = (alpha = mat_alpha, sigma = mat_sigma, mu = mat_mu)

menura_sim_exper_mat_OU(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma, mat_mu, tree1)