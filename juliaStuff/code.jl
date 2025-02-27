using Pkg

Pkg.activate(".")
Pkg.develop(path="/home/simoneb/Desktop/JMMenura/")

# Set to location of JMMenura
##include("../JMMenura/src/JMMenura.jl")

# Import packages
using JMMenura
using Phylo, Distributions, Random, JLD2, LinearAlgebra, PosDefManifold, Distances, StatsBase, GpABC

# Number of traits to simulate
n = 4 # SET

# Open tree - Set file path to tree
tree1 = open(parsenewick, "bigsim.tre")

# Find the root node of the trait
root = getroot(tree1)
root_num = tree1.nodedict[root.name]

# Set trait parameters needed to evolve traits
# Should be a vector of length n
trait_alpha = repeat([1.0], n)
trait_mu = repeat([0.0], n) ##
trait_sigma = repeat([sqrt(2)], n)

# Parameters used for reference simulation
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

@load "./P0.jld2"

# Variables needed for OU matrix model
mat_alpha = 1
mat_sigma = sqrt(2)
mat_mu = copy(P0)

# Create matrix dictionary
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

# Specify what evolution functions are used
# Can set dt the evolution lenght
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.01)

@time ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func
, mat_evol_func, 0.0, trait_mu, P0, true);

ref_data = get_data(ref_sim)
@save "ref_data.jld2" ref_data
@save "tree.jld2" tree1

##SP1 = Distribution{Univariate, Continuous}

##struct HalfNormal <: SP1
##    σ::Float64
##    end

##function Base.rand(rng::AbstractRNG, d::HalfNormal)
##    z = randn(rng)
##    return abs(d.σ * z)
##    end
sigmaPrior = 10.0
## prior = HalfNormal(sigmaPrior)
prior = Truncated(Normal(0, sigmaPrior), 0, Inf)
para = JMMABCAlphaDifferentConstant([prior for _ in 1:n], trait_mu, trait_sigma, prior, mat_mu, mat_sigma, n)

@time thresholds = test_threshold(ref_data, tree1, para, zeros(n), P0, 10, dt = 0.01, each = true) ## changed threshold samples to 10
threshold = sort(thresholds)[4] ## 4th highest threshold

# Save threshold
@save "./threshold.jld2" threshold
s_para_result = Vector{Any}()
@save "./s_para_result.jld2" s_para_result

n_particles = 1000
@load "./threshold.jld2" threshold
@load "ref_data.jld2" ref_data
@time run_result = menura_bayesian(ref_data, tree1, para, zeros(n), P0, threshold,
    n_particles, dt = 0.01,
    max_iter = 25000 * n_particles, each = true, verbose=true);

# Save model results
@load "./s_para_result.jld2" s_para_result
push!(s_para_result, run_result)
@save "s_para_result.jld2" s_para_result

exit()

