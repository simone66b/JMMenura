using .JMMenura
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, Profile, Random, LinearAlgebra
# read in data

Pkg.activate(".")
tree_anole = open(parsenewick, "anoles_data//pruned7.tre")

plot(tree_anole)

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]

function read_cov_mat(species)
    cov = open("anoles_data/"*species*".txt","r") do datafile
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

cov_mats = read_cov_mat.(files)

cov_mean = mean(cov_mats)

trait_data = DataFrame(XLSX.readtable("anoles_data/Adult measurements for divergence.xlsx", "Pmatrix Measurements with outli"))

overall_trait_mean = describe(trait_data[:, 3:11], :mean)[1:9, 2]

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

species_traits = [species_subset(trait_data[:,2:11], x) for x in names]

trait_means = [describe(df[:,2:10], :mean)[2:9, 2] for df in species_traits]

trait_sd = [describe(df[:,2:10], :std)[1:9, 2] for df in species_traits]
overall_trait_sd = mean(trait_sd)

cov_sd = std(cov_mats)

# Set up parameters
prior = Uniform(0,3)
para = JMMABCAlphaEqualConstant(prior, overall_trait_mean[2:end], overall_trait_sd[2:end], prior, cov_mean, cov_sd, 8)

# get_reference_data
data = [reduce(hcat, trait_means)..., reduce(hcat, cov_mats)...]
ref_data = reshape(data, length(data), 1)

##########################
# Theoretical parameters #
##########################

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
trait_parameters = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha = 1 .* Matrix(1I, n, n)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(13 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol(dt = 0.041)
trait_evol_func = trait_evol(dt = 0.041)

##############
# Simulation #
##############

thresholds = test_threshold(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, 2000, dt = 0.041)

# @time thresholds = menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, 1000.0, 2000, dt = 0.041).distances[1]

histogram(thresholds)

threshold = percentile(thresholds, 0.1)

# warm up 
@time menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, threshold, 30, dt = 0.041)

@profview menura_bayesian(ref_data, tree1, para, overall_trait_mean[2:end], cov_mean, threshold, 1)

@profview menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, threshold, 1)

out = menura_bayesian(ref_data, tree1, para, overall_trait_mean[2:end], cov_mean, threshold, 10)


