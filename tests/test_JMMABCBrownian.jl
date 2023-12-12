using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances
pyplot()


Random.seed!(1)

# number of parameters
n = 8

# Creating tree
tree = Ultrametric(7)
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
trait_parameters_true = Dict(13 => (alpha = alpha1, mu = mu1, sigma = sigma1))
trait_parameters = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha = 0 .* Matrix(1I, n, n)
mat_sigma = (1 / sqrt(2) ) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(13 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

para_priors = [Uniform(0,3), Gamma(2,0.25)]

ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

ther_ref_data = get_data(ther_ref_sim[1])

parameters = JMMABCBrownian(Uniform(0,3), mu1, sigma1, mat_mu, Gamma(2, 0.25), n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, mu1, P0, 100)

histogram(thresholds)

@time out = menura_bayesian(ther_ref_data, tree1, parameters, mu1, P0, 245.0, 40)

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, mu1, P0, 245.0, 40)


x = out.population[:,1]
y = out.population[:,2]
append!(x, out2.population[:,1])
append!(y, out2.population[:,2])
histogram2d(x, y, normalize = :pdf, show_empty_bins = true)