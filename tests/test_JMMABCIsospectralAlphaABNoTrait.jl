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
P0 = rand(Wishart(100, Matrix(1I, n, n)  ))

# create trait dictionary
trait_parameters_true = Dict(13 => nothing)
trait_parameters = nothing

# Variables needed for OU matrix model
mat_a = 1
mat_b = 1

# create matrix dictionary
mat_parameters_true = Dict(13 => (a = 1.0,  b = 1.0))
mat_parameters = (a = mat_a, b = mat_b)

mat_evol_func = mat_evol_isospectral()

a_prior = Uniform(1, 10)
b_prior = Uniform(1, 10)

ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, nothing, mat_evol_func, 0.0, [0.0], P0, true)

ther_ref_data = get_data_no_trait(ther_ref_sim)

parameters = JMMABCIsospectralAlphaABNoTrait(a_prior, b_prior, n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, [0.0], P0, 20, each = true, distance_function = mat_distance(4,7))

threshold = percentile(thresholds,5)

@time out = menura_bayesian(ther_ref_data, tree1, parameters, [0.0], P0, threshold, 1, distance_function = mat_distance(4,7))

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, [0.0], P0, threshold, 40, distance_function = mat_distance(4,7))


x = out.population[1][:,1]
y = out.population[1][:,2]
append!(x, out2.population[1][:,1])
append!(y, out2.population[1][:,2])
histogram2d(x, y, normalize = :pdf, show_empty_bins = true)