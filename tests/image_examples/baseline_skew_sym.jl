# Creates a basic evolution plot
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold

pyplot()

Random.seed!(1)

# number of parameters
n = 6

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

# Variables needed for skew symmetric matrix model
mat_a = 1
mat_b = 1

# create matrix dictionary
mat_parameters = Dict(11 => (a = mat_a, b = mat_b))

mat_evol_func = mat_evol_isospectral()
trait_evol_func = trait_evol()

menura_parameter_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, false)

plot_data(tree1, 1, 2, ylim = (-2, 2), zlim = (-2.0, 2.0), legend = false, reuse = false)

for node in getnodes(tree1)
    @show eigen(node.data["mat_trace"][1]).values
end

