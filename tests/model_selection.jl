using .JMMenura
using GpABC, Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold


pyplot()

Random.seed!(1)

# number of parameters
n = 4

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))

# traits needed to evolve traits
trait_alpha = repeat([1.0], n)
trait_mu = repeat([0.0], n)
trait_sigma = repeat([1.0], n)

# create trait dictionary
trait_parameters = Dict(11 => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

# Variables needed for OU matrix model
mat_alpha = 1 .* ones(n,n)
mat_mu = copy(P0)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)

# create matrix dictionary
mat_parameters = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

reference_data = menura_para_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, trait_mu, P0, true)

plot_data(tree1, 1, 2, ylim = (-2, 2), zlim = (-2.0, 2.0), legend = false, reuse = false)

###################
# Model selection #
###################

OU_prior = Gamma(2, 0.25)

OU_para = JMMABCAlphaEqualConstant(OU_prior, trait_mu, trait_sigma, OU_prior, mat_mu, mat_sigma, n)

BW_prior = Normal(0, 0)

BW_para = JMMABCAlphaEqualConstant(BW_prior, trait_mu, trait_sigma, BW_prior, mat_mu, mat_sigma, n)

Iso_prior = Gamma(2, 0.25)

Iso_para = JMMABCIsospectralAlpha(Iso_prior, trait_mu, trait_sigma, 3, 3, n)

OU_func = create_bayesian_sim(tree1, OU_para, trait_mu, mat_mu, dt = 0.01)

BW_func = create_bayesian_sim(tree1, BW_para, trait_mu, mat_mu, dt = 0.01)

Iso_func = create_bayesian_sim(tree1, Iso_para, trait_mu, mat_mu, dt = 0.01)

model_sim_functions = [OU_func, BW_func, Iso_func]

priors = [get_priors(OU_para), get_priors(BW_para), get_priors(Iso_para)]

dist_func = trait_mat_distance(n, 6, err_thres = 10^-14)

ref_data = get_data(reference_data[1])

out = SimulatedModelSelection(ref_data, model_sim_functions, priors, [29.0], 400, distance_function = dist_func,
max_iter = 1000)


#################################
# Isospectral model to generate #
#################################

iso_mat_evol_func = mat_evol_skew_symmetric()

mat_a = 1
mat_b = 1

# create matrix dictionary
iso_mat_parameters = Dict(11 => (a = mat_a, b = mat_b))

iso_reference_data = menura_para_descend!(iso_mat_parameters, trait_parameters, tree1, trait_evol_func, iso_mat_evol_func, 0.0, trait_mu, P0, true)

iso_ref_data = get_data(iso_reference_data[1])

iso_out = SimulatedModelSelection(iso_ref_data, model_sim_functions, priors, [40.0], 100, distance_function = dist_func,
max_iter = 200000)

# histogram(iso_out.smc_outputs[3].population, normalize = :pdf)

# plot!(Gamma(2, 0.25))

# density!(iso_out.smc_outputs[3].population, trim = true)