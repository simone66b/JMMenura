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
P0 =rand(Wishart(100, Matrix(1I, n, n)  ))

# traits needed to evolve traits
alpha1 = 1.5*rand()*ones(n)
alpha2 = 1.5*rand()*ones(n)
mu1 = repeat([0.0], n)
sigma1 = repeat([1.0], n)

# create trait dictionary
trait_parameters_true = Dict(13 => (alpha = alpha1, mu = mu1, sigma = sigma1), 
                                9 => (alpha = alpha2, mu = mu1, sigma = sigma1))
trait_parameters = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha1 = 1
mat_alpha2 = 2
mat_sigma = 0
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(13 => (alpha = mat_alpha1, mu = mat_mu, sigma = mat_sigma)
                            , 9 => (alpha = mat_alpha2, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol_affine(dt = 0.005)
trait_evol_func = trait_evol(dt =0.005)

α_prior = Uniform(0,4)
ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, ones(n), P0, true)

plot_data(tree1, 1, ylim = (-2, 10), zlim = (-2.0, 2.0), legend = false, reuse = false)


ther_ref_data = get_data(ther_ref_sim)

parameters = JMMABCAlphaDifferentEqual(Dict(13 => α_prior, 9 => α_prior), Dict(13 => mu1, 9 => mu1), Dict(13 => sigma1, 9 => sigma1),
                                        Dict(13 => α_prior, 9 => α_prior), Dict(13 => mat_mu), Dict(13 => mat_sigma), n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, mu1, P0, 20, dt = 0.005)

threshold = quantile(thresholds, (0:2)./length(thresholds))[2]

out = menura_bayesian(ther_ref_data, tree1, parameters, ones(n), P0, threshold, 10, dt = 0.005)

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, ones(n), P0, threshold, 50, dt = 0.005)


x = out.population[1][:,1]
y = out.population[1][:,2]
append!(x, out2.population[1][:,1])
append!(y, out2.population[1][:,2])
histogram2d(x, y, normalize = :pdf, show_empty_bins = true)

histogram(x, normalize = :pdf, legend = false)