using Pkg
Pkg.activate(".")
include("/home/simoneb/Desktop/JMMenura/src/JMMenura.jl")

## cd("/home/simoneb/Desktop/JMMenura")
## Pkg.develop(path="/home/simoneb/Desktop/JMMenura")

using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions
using .JMMenura

function read_cov_mat(species)
    cov = open("/home/simoneb/Desktop/JMMenura/anoles_data/"*species*".txt","r") do datafile # Change as needed
        reduce(hcat,[parse.(Float64, split(line)) for line in eachline(datafile)])
    end
    return cov
end

function species_subset(df, name)
    subset(df, :Species => species -> [coalesce(occursin(name,x), false) for x in species])
end

function kernel(distance, sigma=3.0)
    exp.(- distance.^2.0 ./ (2.0 * sigma^2.0))
end

impTraits(x, y) = abs.(x - y) 

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
    setnodedata!(tree_anole, name, "trait_para", [species_mu])
    setnodedata!(tree_anole, name, "mat_para", [species_cov])
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
N= 5000 ################## Change as required
trait_alpha = repeat([0.0], 8)
trait_mu = [mean(getindex.(trait_means, i)) for i in 1:8] #### repeat([0.0], n) ##
### trait_sigma = repeat([sqrt(2)], n)
sigmaPrior = 50.0
mat_mu = copy(P0)
## mat_sigma = sqrt(2)
## priorvec = repeat([Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)], 8)

prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], 9)
root_num = getroot(tree_anole).id
data = [trait_means..., cov_mats...]
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.01, cond_threshold=1.0e7)

root_num = getroot(tree_anole).id

## prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
##priormat = Truncated(Normal(0.0, 0.5), 0.0, Inf)
## priormat = Uniform(0, 2.0)
##priorvec = repeat([prior], 8)## 8 traits and one for the matrix_diff
## push!(priorvec, priormat)

## sigmasAll = [rand.(priorvec) for i in 1:50000]## 8 + 1 draws from prior. 5000 particles
mat_alpha = 0.0
## number of particles
traits = 8
matrixTraits = traits + 1
species = 7
function sim(N)
    res = []
    p = Progress(N, desc="Processing: ")  # Initialize progress meter
    # Change as needed
    
    for j in 1:N ## major loop for 5000 particles     
        # Loop to rerun till stability
        sol_stable = false
        result = -1
        sigmasAll = rand.(priorvec)
        ##reruns = 0
        ## max_reruns = 10
        while !sol_stable ## && reruns <= max_reruns
            sigmasAll = rand.(priorvec)
        mat_parameters_true = Dict(root_num => (alpha = mat_alpha, mu = mat_mu, sigma = sigmasAll[9]))
        trait_parameters_true = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = sigmasAll[1:8]))
        tree_anole1 = open(parsenewick, "/home/simoneb/Desktop/JMMenura/anoles_data/prunedscaled.tre") # Change as needed
            try
            result = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree_anole1, 
            trait_evol_func, mat_evol_func, 0.0, trait_mu, P0, true);
            catch e
            ## reruns += 1
            continue ### sol_stable = false 
            end
         sol_stable = result[2] 
        end
        rundat = get_data2(result)

        dattraits = rundat[1:7] # 7 species
        reftraits = data[1:7] # 7 species
        vals = impTraits.(dattraits, reftraits)
        traitImportances = sum(kernel.(vals))
        datmats = rundat[8:14]
        refmats =data[8:14]
        datmats1, refmats1 = Hermitian.(datmats), Hermitian.(refmats)
        matImportances = sum(kernel.(sqrt.(distanceSqr.(Fisher, datmats1, refmats1))))
        Importances = [sigmasAll, push!(traitImportances, matImportances)]
        push!(res, Importances)
        next!(p)
    end
    res
end

tst = sim(5000)
sigmas = [tst[i][1] for i in 1:N]
wts = [tst[i][2] for i in 1:N]
@save "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/BManoles50.jld2" sigmas wts
# for k in 1:9
# sigmastraitk = [sigmas[i][k] for i in 1:N]
# wtstraitk = [wts[i][k] for i in 1:N]
# wtstrait1normalised = Weights(wtstraitk)

# samps1 = sample(sigmastraitk, wtstrait1normalised, 10000, replace=true)
# his = histogram(samps1, density=true)

# if k == 9 
#     x = range(0,3, length=100)
# else
# x = range(0, 4, length=100)
# end
# plot!(x, pdf.(priorvec[k], x), color=:red)
# end
# show(his)