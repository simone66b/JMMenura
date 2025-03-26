using Pkg
  
Pkg.activate(".")
include("/home/simoneb/Desktop/JMMenura/src/JMMenura.jl")

## cd("/home/simoneb/Desktop/JMMenura")
## Pkg.develop(path="/home/simoneb/Desktop/JMMenura")

using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions

using .JMMenura
#### using Debugger

function read_cov_mat(species)
    cov = open("/home/simoneb/Desktop/JMMenura/anoles_data/"*species*".txt","r") do datafile # Change as needed
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

function kernel(distance, sigma=500)
    exp.(- distance.^2 ./ (2 * sigma^2))
end

impTraits(x, y) = (x - y) .^ 2.0

trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))  # Change as needed


tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed
files = ["cris", "pulc", "grah", "ever", "line", "sagr", "smar"]
names2 = ["cristatellus", "pulchellus", "evermanni", "grahami", "lineatopus", "sagrei", "smaragdinus"]


cov_mats = []
for (i, name) in enumerate(names2)
    species_traits = species_subset(trait_data[:,2:11], name)
    species_mu = describe(species_traits[:,2:10], :mean)[2:9, 2]
    species_cov = read_cov_mat(files[i])
    push!(cov_mats, species_cov)
    setnodedata!(tree_anole, name, "trait_trace", [species_mu])
    setnodedata!(tree_anole, name, "mat_trace", [species_cov])
end

species_traits = [species_subset(trait_data[:,2:11], x) for x in names2]

trait_means = [describe(df[:,2:10], :mean)[2:9, 2] for df in species_traits]
## cov_mats = read_cov_mat.(files)

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

## trait_alpha = repeat([1.0], n)
trait_mu = trait_means[1]  ##repeat([0.0], n) ##
trait_alpha = repeat([0.0], n)
sigmaPrior = 10.0
## mat_alpha = 0.5
mat_mu = copy(P0)
## mat_sigma = sqrt(2)
##priors = repeat([Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)], 9)
##root_num = getroot(tree_anole).id
a_sim_res = []
data = [trait_means..., cov_mats...]
##trait_evol_func = trait_evol(dt = 0.01)
##mat_evol_func = mat_evol_affine(dt = 0.01)

trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.01)   


root_num = getroot(tree_anole).id

prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], 9)## 8 traits and one for the matrix_diff

## alphasAll = [rand.(priorvec) for i in 1:5000]## 8 + 1 draws from prior. 5000 particles
sigmasAll = [rand.(priorvec) for i in 1:5000] 

N= 5000 ## number of particles
traits = 8
matrixTraits = traits + 1
species = 7

######################################################################################
function sim(N)
res = []
p = Progress(N, desc="Processing: ")  # Initialize progress meter

j = 1 ##for j in 1:N ## major loop for 5000 particles
while j <= N
try
mat_parameters_true = Dict(root_num => (alpha = 0.0, mu = mat_mu, sigma = sigmasAll[j][9]))
trait_parameters_true = Dict(root_num => (alpha = repeat([0.0], 8), mu = trait_mu, sigma = sigmasAll[j][1:8]))

tree_anole = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed
    
result = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree_anole, trait_evol_func, 
mat_evol_func, 0.0, trait_mu, P0, true);

rundat = get_data2(result)

dattraits = rundat[1:7] # 7 species
reftraits = data[1:7] # 7 species

vals = impTraits.(dattraits, reftraits)
traitImportances = sum(kernel.(vals))

datmats = rundat[8:14]
refmats =data[8:14]
datmats1, refmats1 = Hermitian.(datmats), Hermitian.(refmats)
matImportances = sum(kernel.(distanceSqr.(Fisher, datmats1, refmats1)))

Importances = [sigmasAll[j], push!(traitImportances, matImportances)]
push!(res, Importances)
println("Finished iteration $j...")
j += 1
## next!(p)
catch error
end

end
return res
end

tst = sim(500)
sigmas = [tst[i][1] for i in 1:79]
wts = [tst[i][2] for i in 1:79]

## for k in 1:9
# alphastraitk = [alphas[i][k] for i in 1:5000]
# wtstraitk = [wts[i][k] for i in 1:5000]
# wtstrait1normalised = Weights(wtstraitk)

# samps1 = sample(alphastraitk, wtstrait1normalised, 100000, replace=true)
# his = histogram(samps1, normalize=true)

# x = range(0, 40, length=100)

# plot!(x, pdf.(prior, x), color=:red)
# display(his)
## end