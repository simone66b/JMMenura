using .JMMenura
using Phylo, Distributions, Pkg, Plots
# read in data

Pkg.activate(".")
tree = open(parsenewick, "anoles_data//pruned7.tre")

plot(tree)

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]

function read_cov_mat(species)
    cov = open("anoles_data/"*species*".txt","r") do datafile
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

cov_mats = read_cov_mat.(files)

cov_mean = mean(cov_mats)