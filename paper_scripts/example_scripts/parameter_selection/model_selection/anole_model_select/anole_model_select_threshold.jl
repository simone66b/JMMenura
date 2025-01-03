using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, GpABC
Pkg.activate(".")

include("./../../../JMMenura/src/JMMenura.jl")

using .JMMenura

tree_anole = open(parsenewick, "./../..//anoles_data//prunedscaled.tre")

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]

function read_cov_mat(species)
    cov = open("./../..//anoles_data/"*species*".txt","r") do datafile
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

cov_mats = read_cov_mat.(reverse(files))

cov_mean = mean(cov_mats)

trait_data = DataFrame(XLSX.readtable("./../..//anoles_data/Adult measurements for divergence.xlsx", "Pmatrix Measurements with outli"))

overall_trait_mean = describe(trait_data[:, 3:11], :mean)[1:9, 2]

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

species_traits = [species_subset(trait_data[:,2:11], x) for x in reverse(names)]

trait_means = [describe(df[:,2:10], :mean)[2:9, 2] for df in species_traits]

trait_sd = [describe(df[:,2:10], :std)[1:9, 2] for df in species_traits]
overall_trait_sd = mean(trait_sd)

cov_sd = std(cov_mats)


n = 8

trait_sigma = repeat([sqrt(2)], n)

mat_sigma = sqrt(2)

#####################
# Set up parameters #
#####################

OU_prior = Gamma(2, 0.25)

OU_para = JMMABCAlphaEqualConstant(OU_prior, overall_trait_mean[2:end], trait_sigma, OU_prior, cov_mean, mat_sigma, 8)

BW_prior = Gamma(2, 0.25)

BW_para = JMMABCBrownian(BW_prior, overall_trait_mean[2:end], trait_sigma, cov_mean, BW_prior, 8)

Iso_prior = Gamma(2, 0.25)

Iso_para = JMMABCIsospectralAlphaAB(Iso_prior, overall_trait_mean[2:end], trait_sigma, Iso_prior, Iso_prior, 8)



OU_func = create_bayesian_sim(tree_anole, OU_para, overall_trait_mean[2:end], cov_mean, dt = 0.005, each = true)

BW_func = create_bayesian_sim(tree_anole, BW_para, overall_trait_mean[2:end], cov_mean, dt = 0.005, each = true)

Iso_func = create_bayesian_sim(tree_anole, Iso_para, overall_trait_mean[2:end], cov_mean, dt = 0.005, each = true)

model_sim_functions = [OU_func, BW_func, Iso_func]

priors = [get_priors(OU_para), get_priors(BW_para), get_priors(Iso_para)]

dist_func = trait_mat_distance(8, 7)
######################
# get_reference_data #
######################

data = [reduce(hcat, trait_means)..., reduce(hcat, cov_mats)...]
ref_data = reshape(data, length(data), 1)

######################
# Perform simulation #
######################

n_threshold = 1000

thresholds = [test_threshold(ref_data, tree_anole, OU_para, overall_trait_mean[2:end], cov_mean, n_threshold, dt = 0.005, each = true), 
                test_threshold(ref_data, tree_anole, BW_para, overall_trait_mean[2:end], cov_mean, n_threshold, dt = 0.005, each = true), 
                test_threshold(ref_data, tree_anole, Iso_para, overall_trait_mean[2:end], cov_mean, n_threshold, dt = 0.005, each = true)]

threshold = minimum([sort(model_thresholds)[4] for model_thresholds in thresholds])

@save "./threshold.jld2" threshold

a_model_result = Vector{Any}()

@save "./a_model_result.jld2" a_model_result
