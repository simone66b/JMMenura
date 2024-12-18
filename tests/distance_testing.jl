# Creates a basic evolution plot
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances


# testing how the distance functions work

# P0 = rand(Wishart(100, Matrix(1I, n, n)  ))
# P1 = rand(Wishart(100, Matrix(1I, n, n)  ))


# P0H = Hermitian(P0)
# P1H = Hermitian(P1)

# eigen(P0H)
# eigen(P1H)

# logMap(Fisher, P0H, P1H)
# PosDefManifold.distance(Fisher, P0H, P1H)

# t1 = rand(9)
# t2 = rand(9)

# euclidean(t1, t2)

Random.seed!(1)

# number of parameters
n = 7

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))

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

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

#Approximating 5% cutoff

points = zeros(300)

for i in 1:100
menura_para_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

data1 = get_data(tree1)

menura_para_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

data2 = get_data(tree1)

points[i] = trait_mat_distance(7, 6)(data1, data2)

end

#

histogram(points)