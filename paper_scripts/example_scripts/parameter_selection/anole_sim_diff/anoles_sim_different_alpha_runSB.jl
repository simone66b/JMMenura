using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
## Pkg.activate(".")

cd("/home/simoneb/Desktop/JMMenura/")

## include("/home/simoneb/Desktop/JMMenura/src/JMMenura.jl")

using JMMenura

tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]

function read_cov_mat(species)
    cov = open("/home/simoneb/Desktop/JMMenura/anoles_data/"*species*".txt","r") do datafile # Change as needed
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))  # Change as needed

for (i, name) in enumerate(names)
    species_traits = species_subset(trait_data[:,2:11], name)
    species_mu = describe(species_traits[:,2:10], :mean)[2:9, 2]

    species_cov = read_cov_mat(files[i])
    setnodedata!(tree_anole, name, "trait_trace", [species_mu])
    setnodedata!(tree_anole, name, "mat_trace", [species_cov])
end

species_traits = [species_subset(trait_data[:,2:11], x) for x in reverse(names)]

trait_means = [describe(df[:,2:10], :mean)[2:9, 2] for df in species_traits]


n = 8 ## 8 traits
GancVec = [0.285, 0.118, 0.131, 0.053, 0.212, 0.172, 0.188, 0.210, 0.277,
            0.230, 0.103, 0.074, 0.073, 0.045, 0.023, 0.800, 0.317, 0.054,
            0.121, 0.152, 0.147, 0.654, 0.028, 0.047, 0.124, 0.148, 0.858,
            0.749, 0.581, 0.626, 0.906, 0.553, 0.644, 0.708, 0.731, 0.872]
Ganc =zeros(8, 8)
idx = tril!(trues(size(Ganc)))
Ganc[idx] .= GancVec
Ganc .= Hermitian(Ganc, :L)
P0 = Ganc #### FullData[7, :Gmatrix]; ## randPosDefMat(8)

trait_alpha = repeat([1.0], n)
trait_mu = repeat([0.0], n) ##
trait_sigma = repeat([sqrt(2)], n)
sigmaPrior = 10.0
mat_alpha = 0.5
## @load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"
mat_mu = copy(P0)
mat_sigma = sqrt(2)


root_num = getroot(tree_anole).id
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)

data = [trait_means..., cov_mats...]
## ref_data = reshape(data, length(data), 1)
##trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_means, sigma = trait_sigma))
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.0)

## mat_mu = copy(P0)
tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed

@time result = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree_anole, trait_evol_func, 
mat_evol_func, 0.0, trait_mu, P0, true);

rundat = get_data2(result)

dat = rundat
refdat = data
function importance(dat, refdat)
end

    dattraits = dat[1:7] # 7 species
    reftraits = refdat[1:7] # 7 species

test(x, y) = 1.0 ./ (x - y) .^ 2.0
vals = test.(dattraits, reftraits)
traitImportances = reduce(.+, vals)

m = SymmetricPositiveDefinite(8)
testmat(x, y) = 1.0 ./ fisher_rao_distance(x, y) .^ 2.0
datmats = dat[8:14]
refmats = refdat[8:14]
testmat.(datmats, refmats)
testmat(x, y) = 1.0 ./ fisher_rao_distance(x, y) .^ 2.0

@load "./a_sim_diff_result.jld2" a_sim_diff_result

push!(a_sim_diff_result, run_result)

@save "a_sim_diff_result.jld2" a_sim_diff_result

