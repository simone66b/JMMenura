using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2
Pkg.activate(".")

include("./../../../../../JMMenura/src/JMMenura.jl")

using .JMMenura

tree_anole = open(parsenewick, "./../../../..//anoles_data//prunedscaled.tre")

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]

function read_cov_mat(species)
    cov = open("./../../../..//anoles_data/"*species*".txt","r") do datafile
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

cov_mats = read_cov_mat.(reverse(files))

cov_mean = mean(cov_mats)

trait_data = DataFrame(XLSX.readtable("./../../../..//anoles_data/Adult measurements for divergence.xlsx", "Pmatrix Measurements with outli"))


overall_trait_mean = describe(trait_data[:, 3:11], :mean)[1:9, 2]

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

species_traits = [species_subset(trait_data[:,2:11], x) for x in reverse(names)]

trait_means = [describe(df[:,2:10], :mean)[2:9, 2] for df in species_traits]

trait_sd = [describe(df[:,2:10], :std)[1:9, 2] for df in species_traits]
overall_trait_sd = mean(trait_sd)

cov_sd = std(cov_mats)

#####################
# Set up parameters #
#####################

trait_sigma = repeat([sqrt(2)], 8)
mat_sigma = sqrt(2)

prior = Gamma(2, 0.25)
para = JMMABCAlphaDifferentConstant([prior for _ in 1:8], overall_trait_mean[2:end], trait_sigma, prior, cov_mean, mat_sigma, 8)

######################
# get_reference_data #
######################

data = [reduce(hcat, trait_means)..., reduce(hcat, cov_mats)...]
ref_data = reshape(data, length(data), 1)

###############
# Warm up sim #
###############

@load "./threshold.jld2" threshold

n_particles = 200

run_result = menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, threshold, n_particles, dt = 0.005, max_iter = 50000*n_particles, each = true)

@load "./a_sim_diff_result.jld2" a_sim_diff_result

push!(a_sim_diff_result, run_result)

@save "a_sim_diff_result.jld2" a_sim_diff_result

