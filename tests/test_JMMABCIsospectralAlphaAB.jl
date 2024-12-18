using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances, StatsBase
pyplot()


Random.seed!(1)

# number of parameters
n = 4

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
mat_a = 1
mat_b = 1

# create matrix dictionary
mat_parameters_true = Dict(13 => (a = 1.0,  b = 1.0))
mat_parameters = (a = mat_a, b = mat_b)

mat_evol_func = mat_evol_isospectral()
trait_evol_func = trait_evol()

α_prior = Uniform(0,3)
a_prior = Uniform(1, 10)
b_prior = Uniform(1, 10)

ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

ther_ref_data = get_data(ther_ref_sim[1])

parameters = JMMABCIsospectralAlphaAB(α_prior, mu1, sigma1, a_prior, b_prior, n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, mu1, P0, 2000, each = false)

threshold = percentile(thresholds,0.1)

@time out = menura_bayesian(ther_ref_data, tree1, parameters, mu1, P0, threshold, 100)

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, mu1, P0, 245.0, 40)


x = out.population[:,1]
y = out.population[:,2]
append!(x, out2.population[:,1])
append!(y, out2.population[:,2])
histogram2d(x, y, normalize = :pdf, show_empty_bins = true)