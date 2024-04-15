using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances
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
mat_alpha = 1 .* Matrix(1I, n, n)
mat_sigma = 1 / sqrt(2) 
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(13 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol_affine(dt = 0.05)

ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, nothing, mat_evol_func, 0.0, [0.0], P0, true)

ther_ref_data = get_data_no_trait(ther_ref_sim)

parameters = JMMABCAlphaEqualConstantNoTrait(Uniform(0,3), mat_mu, mat_sigma, n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, [0.0], P0, 10, distance_function = mat_distance(4,7), each = true)

histogram(thresholds)

threshold = sort(thresholds)[2]

@time out = menura_bayesian(ther_ref_data, tree1, parameters, [0.0], P0, threshold, 1, distance_function = mat_distance(4,7), each = true)

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, [0.0], P0, threshold, 40, distance_function = mat_distance(4,7), each = true)


x = out.population[1][:,1]
append!(x, out2.population[1][:,1])
histogram(x, normalize = :pdf)