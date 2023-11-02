using .JMMenura
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2

###########################
# Read in and format data #
###########################

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

#####################
# Set up parameters #
#####################

prior = Gamma(2, 0.25)
para = JMMABCAlphaEqualConstant(prior, overall_trait_mean[2:end], overall_trait_sd[2:end], prior, cov_mean, cov_sd, 8)

######################
# get_reference_data #
######################

data = [reduce(hcat, trait_means)..., reduce(hcat, cov_mats)...]
ref_data = reshape(data, length(data), 1)

######################
# Perform simulation #
######################

thresholds = test_threshold(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, 100)

histogram(thresholds)

threshold = percentile(thresholds, 2)

@time out = menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, [threshold], 10, dt = 0.05)

out2 = menura_bayesian(ref_data, tree_anole, para, overall_trait_mean[2:end], cov_mean, [threshold], 1000, dt = 0.05)


pyplot()

x = out2.population[1][:,1]
y = out2.population[1][:,2]
# append!(x, out2.population[:,1])
# append!(y, out2.population[:,2])
histogram2d(x, y, normalize = :pdf, show_empty_bins = true, xlab = "Trait alpha", ylab = "G matrix alpha")
