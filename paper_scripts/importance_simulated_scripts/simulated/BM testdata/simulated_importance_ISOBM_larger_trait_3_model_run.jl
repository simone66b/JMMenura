using Pkg
## cd("/Users/coope/OneDrive/Documents/Uni/Phylogenetics_coding/importance_scripts")
include("/home/simoneb/Desktop/JMMenura/src/JMMenura.jl")
cd("/home/simoneb/Desktop/JMMenura")
Pkg.develop(path="/home/simoneb/Desktop/JMMenura")
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions
using .JMMenura

##################################
# Load reference simulation data #
##################################
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/Test Data sets and trees/BM_larger_trait_3_para_ref_data5D.jld2" para_ref_data_tree

trees = [para_ref_data_tree.ref.mem[i][1] for i in 1:5]
##     open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim$i.tre")
data = [para_ref_data_tree.ref.mem[i][2] for i in 1:5]
#####################
# Set up parameters #
#####################
n = 4 ## traits
num_species = 50 
## tree1 = open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim.tre")
## nu = Ultrametric(num_species)
## trees = rand(nu,  ['A':'E';])
## Phylo.write("/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/fiveTrees.tre", trees)

## trees = open(parsenexus, "fiveTrees.tre")

##trees = [para_ref_data_tree.tre for i in 1:5]
##    open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim$i.tre")

# BM simulation for traits
trait_alpha = repeat([0.0], n)
trait_mu = repeat([0.0], n) ##

# BM simulation for matrix
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/P1.jld2"

# mat_alpha = 0.0
# mat_a = 5
# mat_b = 5
# mat_mu = copy(P0)

# data = para_ref_data
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_isospectral(dt = 0.01)

# Same starting conditions as OU
start_trait_alpha = rep([0.0], n) ## [2, 4, 6, 8]
start_trait_mu = repeat([0.0], n)
start_trait_sigma = repeat([sqrt(2)], n)

trait_start = repeat([0.0], n) ## start_trait_mu + 3*(start_trait_sigma ./ sqrt.(2*start_trait_alpha))
mat_start = P1
mat_mu = copy(P0)
mat_alpha = 0.0
mat_sigma = sqrt(2.0)

sigmaPrior = 50
abPrior = 50
sigmaprior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
sigmapriorvec = repeat([sigmaprior], n) # 4 traits
sigmasAll = [rand.(sigmapriorvec) for i in 1:5000]## 8 + 2 draws from prior. 5000 particles

prior = Truncated(Normal(0.0, abPrior), 0.0, Inf)
priorvecab = repeat([prior], 2) # two for the matrix_diff
abAll = [rand.(priorvecab) for i in 1:5000]## 8 + 2 draws from prior. 5000 particles



function get_data2(sim_data)
    tree = sim_data[1]
    traits = [tip.data["trait_trace"][end] for tip in getleaves(tree)]
    mats = [tip.data["mat_trace"][end] for tip in getleaves(tree)]
    data = [traits..., mats...]
    return data
end

function kernel(distance, sigma=3)
    exp.(- distance.^2 ./ (2 * sigma^2))
end

function impTraits(x, y) 
    return abs.(x - y)
end


function sim(N, tree, data)
# tree= trees[1]
# data = data[1]
# N = 1
    root_num = getroot(tree).id
    res = []
    unrun = []
    p = Progress(N, desc="Processing: ")  # Initialize progress meter
    max_reruns = 500000 # Change as needed
    for j in 1:N ## major loop for 5000 particles
        # Loop to rerun till stability
     ### j = 1
        reruns = 0
        sol_stable = false
        rerun = false
        result = nothing
       while !sol_stable && reruns < max_reruns
            if rerun
                push!(unrun, (sigmasAll[j], abAll[j]))
                sigmasAll[j] = rand.(sigmapriorvec)
                abAll[j] = rand.(priorvecab)
            end
            mat_parameters = Dict(root_num => (a = abAll[j][1], b = abAll[j][2]))
            trait_parameters = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma = sigmasAll[j][1:n]))
            try
            result = menura_parameter_descend!(mat_parameters, trait_parameters, tree, trait_evol_func, mat_evol_func, 0.0, trait_start, mat_start, true)
            sol_stable = true
            catch
                rerun = true
                reruns += 1
                sol_stable = false
            end
            if reruns >= max_reruns
               throw("Number of reruns $reruns exceeding maximum $max_reruns. Try increasing max_reruns.")
            end
        end
        rundat = get_data2(result)
        dattraits = rundat[1:num_species] # 50 species
        reftraits = data[1:num_species] # 50 species
        vals = impTraits.(dattraits, reftraits)
        traitImportances = sum(kernel.(vals))
        datmats = rundat[(num_species+1):(2*num_species)]
        refmats =data[(num_species+1):(2*num_species)]
        datmats1, refmats1 = Hermitian.(datmats), Hermitian.(refmats)
        matImportances = sum(kernel.(sqrt.(distanceSqr.(Fisher, datmats1, refmats1))))
        Importances = push!(traitImportances, matImportances)
        push!(res, Importances)
        next!(p)
    end
    res, unrun
end

# Perform simulation
for i in 1:5 
    tree = trees[i]
    this_data = data[i]
    println("Running simulation for tree $i")
tst, unrun = sim(5000, tree, this_data)
@save "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/Test Data sets and trees/ISOmodelled$i.BMtestdata.jld2" tst unrun
    println("Simulation for tree $i completed and saved.")
end
