using Pkg

Pkg.activate(".")

include("./../../../JMMenura/src/JMMenura.jl")
using .JMMenura
using Phylo, Distributions, Random, JLD2, LinearAlgebra, PosDefManifold, Distances, StatsBase, GpABC

#####################
# Set up parameters #
#####################

n = 4

# Creating tree
tree1 = open(parsenewick, "./../..//anoles_data//bigsim.tre")

time_tot = 1.0
tspan = (0.0, time_tot)

# Get root number
root = getroot(tree1)
root_num = tree1.nodedict[root.name]

# G matrix
@load "./../../anoles_data/P0.jld2"

# traits needed to evolve traits
trait_alpha = repeat([1.0], n)
trait_mu = repeat([0.0], n)
trait_sigma = repeat([sqrt(2)], n)

# create trait dictionary
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))
trait_parameters = (mu = trait_mu, sigma = trait_sigma)

############################
# Set up OU reference data #
############################

# Variables needed for OU matrix model
mat_alpha = 1
mat_sigma = sqrt(2)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol_affine(dt = 0.005)
trait_evol_func = trait_evol(dt = 0.005)

OU_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, trait_mu, P0, true)

OU_ref_data = get_data(OU_ref_sim)

@save "big_OU_ref_data.jld2" OU_ref_data

@save "tree.jld2" tree1


#################
# Create Models #
#################

OU_prior = Gamma(2, 0.25)

OU_para = JMMABCAlphaEqualConstant(OU_prior, trait_mu, trait_sigma, OU_prior, mat_mu, mat_sigma, n)

BW_prior = Gamma(2, 0.25)

BW_para = JMMABCBrownian(BW_prior, trait_mu, trait_sigma, mat_mu, BW_prior, n)

Iso_prior = Gamma(2, 0.25)

Iso_para = JMMABCIsospectralAlphaAB(Iso_prior, trait_mu, trait_sigma, Iso_prior, Iso_prior, n)



OU_func = create_bayesian_sim(tree1, OU_para, trait_mu, mat_mu, dt = 0.005, each = true)

BW_func = create_bayesian_sim(tree1, BW_para, trait_mu, mat_mu, dt = 0.005, each = true)

Iso_func = create_bayesian_sim(tree1, Iso_para, trait_mu, mat_mu, dt = 0.005, each = true)

model_sim_functions = [OU_func, BW_func, Iso_func]

priors = [get_priors(OU_para), get_priors(BW_para), get_priors(Iso_para)]

dist_func = trait_mat_distance(n, 50)

###########################
# Perform Model Selection #
###########################

n_threshold = 1000

OU_thresholds = [test_threshold(OU_ref_data, tree1, OU_para, trait_mu, mat_mu, n_threshold, each = true),
                test_threshold(OU_ref_data, tree1, BW_para, trait_mu, mat_mu, n_threshold, each = true), 
                test_threshold(OU_ref_data, tree1, Iso_para, trait_mu, mat_mu, n_threshold, each = true)]

threshold = minimum([sort(model_thresholds)[4] for model_thresholds in OU_thresholds])

@save "./threshold.jld2" threshold

s_model_big_OU_result = Vector{Any}()

@save "./s_model_big_OU_result.jld2" s_model_big_OU_result
