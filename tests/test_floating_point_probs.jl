# Test file to test floating point errors
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

# number of parameters
n = 9

# Set up trees and time span
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
Random.seed!(1)
tree2 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
Random.seed!(4)
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) # Change to Wishart Distribution

# Looking for small negative eigenvalue
eigen(P0)

# trait vairables - note alpha*mu >= sigma
alpha1 = repeat([1.0], n)
mu1 = repeat([1.0], n) ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n)
parms= (alpha=alpha1, sigma=sigma1)

# Matrix variables
mat_alpha = 0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)

# This should error by an inexact error
data1 = menura_sim_mat_OU_each(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma, mat_mu, tree1)

# Adding a small error value
ϵ = Matrix((10^-15)I, n, n)
P1 = P0 + ϵ
mat_mu1 = copy(P1)

# Check eigenvalues
eigen(P1)

# Should not error now
data2 = menura_sim_mat_OU_each(alpha1, sigma1, mu1, P1, mat_alpha, mat_sigma, mat_mu1, tree1)

plot_data(data2[1], 1, 2)