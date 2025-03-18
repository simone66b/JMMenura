using Pkg
cd("/home/simoneb/Desktop/JMMenura/")
Pkg.activate(".")
## Pkg.develop(path="/home/simoneb/Desktop/JMMenura/")
using Phylo, Distributions, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter
using JMMenura

function read_cov_mat(species)
    cov = open("/home/simoneb/Desktop/JMMenura/anoles_data/"*species*".txt","r") do datafile # Change as needed
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

function kernel(x, sigma=0.5)
    exp(-0.5 * (x/sigma)^2.0)
end


tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed

files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]



trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))  # Change as needed

for (i, name) in enumerate(names)
    species_traits = species_subset(trait_data[:,2:11], name)
    species_mu = describe(species_traits[:,2:10], :mean)[2:9, 2]
    species_cov = read_cov_mat(files[i])
    setnodedata!(tree_anole, name, "trait_trace", [species_mu])
    setnodedata!(tree_anole, name, "mat_trace", [species_cov])
end

## species_traits = [species_subset(trait_data[:,2:11], x) for x in reverse(names)]
## trait_means = [describe(df[:,2:9], :mean)[2:9, 2] for df in species_traits]
cov_mats = read_cov_mat.(files)

n = 8 ## 8 traits
GancVec = [0.285, 0.118, 0.131, 0.053, 0.212, 0.172, 0.188, 0.210, 0.277,
            0.230, 0.103, 0.074, 0.073, 0.045, 0.023, 0.800, 0.317, 0.054,
            0.121, 0.152, 0.147, 0.654, 0.028, 0.047, 0.124, 0.148, 0.858, 
            0.749, 0.581, 0.626, 0.906, 0.553, 0.644, 0.708, 0.731, 0.872]
Ganc =zeros(8, 8)
idx = tril!(trues(size(Ganc)))
Ganc[idx] = GancVec
Ganc = Hermitian(Ganc, :L)
P0 = Matrix(Ganc)

##trait_alpha = repeat([1.0], n)
trait_mu = repeat([0.0], n) ##
trait_sigma = repeat([sqrt(2)], n)
sigmaPrior = 10.0
## mat_alpha = 0.5
## @load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"
mat_mu = copy(P0)
mat_sigma = sqrt(2)
priors = repeat([Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)], 9)
root_num = getroot(tree_anole).id
a_sim_res = []
data = [trait_means..., cov_mats...]
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.01)

p = Progress(5000, desc="Processing: "); 

########################################################################################
for i in 1:5000

alphas = rand.(priors)
mat_alpha  = alphas[9]
trait_alpha = alphas[1:8]
mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = trait_sigma))

## ref_data = reshape(data, length(data), 1)
## trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_means, sigma = trait_sigma))

tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed

result = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree_anole, trait_evol_func, 
mat_evol_func, 0.0, trait_mu, mat_mu, true);

rundat = get_data2(result)

dattraits = rundat[1:7] # 7 species
reftraits = data[1:7] # 7 species
traitDistances = dattraits - reftraits
absTraitDistances = [[abs(x) for x in sub] for sub in traitDistances]
traitWeights = [[kernel(x, 10) for x in sub] for sub in absTraitDistances]
traitWeights2 = reduce(.+, traitWeights)
####function matImportance(datmat, refmat)
    ###abs(PosDefManifold.distance(Fisher, Hermitian(datmat))) / 
   ### abs(PosDefManifold.distance(Fisher, Hermitian(refmat)))
### end

<<<<<<< Updated upstream
vals = impTraits.(dattraits, reftraits)
traitImportances = sum(kernel.(vals))

datmats = rundat[8:14]
refmats =data[8:14]
datmats1, refmats1 = Hermitian.(datmats), Hermitian.(refmats)
matImportances = sum(kernel.(distanceSqr.(Fisher, datmats1, refmats1)))

Importances = [alphasAll[j], push!(traitImportances, matImportances)]
push!(res, Importances)
next!(p)
end
res
end

tst = sim(5000)

=======
datmats = dat[8:14]
refmats = refdat[8:14]
matDistances = abs.(PosDefManifold.distance.(Fisher, Hermitian.(datmats), Hermitian.(refmats)))
matWeights = kernel.(matDistances, 10)
push!(a_sim_res, (alpha=alphas, importance= [traitWeights2..., matWeights...]))
### @load "./a_sim_diff_result.jld2" a_sim_diff_resul
### @save "a_sim_diff_result.jld2" a_sim_diff_result
next!(p)
end
>>>>>>> Stashed changes
