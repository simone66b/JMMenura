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
# eigen(P0) # Checking for small negative eigenvalues

# plot_labelled(tree1) # get number of root node

# traits needed to evolve traits
alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n) ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n)

# create trait dictionary
trait_parameters = Dict(11 => (alpha = alpha1, mu = mu1, sigma = sigma1))

# Variables needed for OU matrix model
mat_alpha = 0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

menura_para_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, mat_mu, false)


