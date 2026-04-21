using Pkg
## cd("/Users/coope/OneDrive/Documents/Uni/Phylogenetics_coding/importance_scripts")
##include("D:/summer project/JMMenura/src/JMMenura.jl")
cd("C:/Users/suqt0/OneDrive/Desktop/UQ/JMMenura")
Pkg.activate(".")
Pkg.instantiate()
## Pkg.develop(path="D:/summer project/JMMenura")
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random
using JMMenura
using Base.Threads

##################################
# Load reference simulation data #
##################################
@load "C:/Users/suqt0/OneDrive/Desktop/UQ/JMMenura/paper_scripts/importance_simulated_scripts/simulated/BM testdata/BM_larger_trait_3_para_ref_data5D.jld2" para_ref_data_tree

trees = [para_ref_data_tree.ref.mem[i][1] for i in 1:5]
##     open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim$i.tre")
data = [para_ref_data_tree.ref.mem[i][2] for i in 1:5]
#####################
# Set up parameters #
#####################
n = 4 ## traits
num_species = 50 
##tree1 = open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim.tre")
## nu = Ultrametric(num_species)
## trees = rand(nu,  ['A':'E';])
## Phylo.write("/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/fiveTrees.tre", trees)

## trees = open(parsenexus, "fiveTrees.tre")

##trees = [para_ref_data_tree.tre for i in 1:5]
##    open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim$i.tre")

# BM simulation for traits
trait_alpha = repeat([0.0], n)
trait_mu = repeat([0.0], n) ##
## trait_sigma = repeat([sqrt(2)], n)

# BM simulation for matrix
@load "C:/Users/suqt0/OneDrive/Desktop/UQ/JMMenura/paper_scripts/example_scripts/anoles_data/P0.jld2"
@load "C:/Users/suqt0/OneDrive/Desktop/UQ/JMMenura/paper_scripts/example_scripts/anoles_data/P1.jld2"
mat_alpha = 0.0
mat_mu = copy(P0)
## mat_alpha = 1.0 
# data = para_ref_data
trait_evol_func = trait_evol(dt = 0.01)
mat_evol_func = mat_evol_affine(dt = 0.01)

# Same starting conditions as OU
start_trait_alpha = repeat([0.0], n) ## [2, 4, 6, 8]
start_trait_mu = repeat([0.0], n)
start_trait_sigma = repeat([sqrt(2)], n)

trait_start = repeat([0.0], n) ## start_trait_mu + 3*(start_trait_sigma ./ sqrt.(2*start_trait_alpha))
mat_start = P1

sigmaPrior = 50  #### cHANGE THIS VALUE TO ADJUST PRIOR
##sigmaPrior = 10 #final optimization
## sigmaPrior = 5
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], n+1) # 4 traits and one for the matrix_diff
sigmasAll = [rand.(priorvec) for i in 1:5000]## 4 + 1 draws from prior. 5000 particles ##
##sigmasAll = [rand.(priorvec) for i in 1:20000]#final optimiztion

##sigmaTrait = 2.5 #final optimization
##sigmaMat   = 5 #final optimization
sigmaTrait = 3.0
sigmaMat = 3.0

function get_data2(sim_data)
    tree = sim_data[1]
    traits = [tip.data["trait_trace"][end] for tip in getleaves(tree)]
    mats = [tip.data["mat_trace"][end] for tip in getleaves(tree)]
    data = [traits..., mats...]
    return data
end

function kernel(distance, sigma=3.0) #final optimization
#function kernel(distance, sigma=3) ## CHANGE SIGMA TO ADJUST KERNEL WIDTH. ALSO CHANGE KERNEL FUNCTION
    exp.(- distance.^2 ./ (2 * sigma^2))
end

function impTraits(x, y) 
    return abs.(x - y)
end

function sim(N, tree, data)

    root_num = getroot(tree).id

    res = []
    unrun = []
    p = Progress(N, desc="Processing: ")  # Initialize progress meter
    max_reruns = 500000 # Change as needed
    reruns = 0
    for j in 1:N ## major loop for 5000 particles
        ##tree1 = open(parsenewick, "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/anoles_data/bigsim.tre")
        # Loop to rerun till stability
        sol_stable = false
        result = -1
        rerun = false
        while !sol_stable && reruns < max_reruns
            if rerun
                push!(unrun, sigmasAll[j])
                sigmasAll[j] = rand.(priorvec)
            end
	    mat_parameters = Dict(root_num => (sigma = sigmasAll[j][n+1], mu = mat_mu, alpha = mat_alpha))
            trait_parameters = Dict(root_num => (alpha = trait_alpha, mu = trait_mu, sigma =sigmasAll[j][1:n]))

            try 
            result = menura_parameter_descend!(mat_parameters, trait_parameters, tree, trait_evol_func, 
            mat_evol_func, 0.0, trait_start, mat_start, true);
            catch
                rerun = true
                reruns += 1
            end
            sol_stable = result[2]
            rerun = true # Only resample if at least one rerun
        end

        if reruns >= max_reruns
            throw("Number of reruns $reruns exceeding maximum $max_reruns. Try increasing max_reruns.")
        end

        rundat = get_data2(result)
        dattraits = rundat[1:num_species] # 50 species
        reftraits = data[1:num_species] # 50 species
        vals = impTraits.(dattraits, reftraits)
        ##traitImportances = mean(kernel.(vals))
        ## traitImportances = mean(kernel.(vals, sigmaTrait))#final optimization
        ## traitImportances = mean(vcat(kernel.(vals)...))
        ##traitImportances = sum(kernel.(vals))
        traitImportances = sum(kernel.(vals, sigmaTrait))
        datmats = rundat[(num_species+1):(2*num_species)]
        refmats =data[(num_species+1):(2*num_species)]
        datmats1, refmats1 = Hermitian.(datmats), Hermitian.(refmats)
        #matImportances = mean(kernel.(sqrt.(distanceSqr.(logEuclidean, datmats1, refmats1)), sigmaMat))#final optimization
        ##matImportances = sum(kernel.(sqrt.(distanceSqr.(Fisher, datmats1, refmats1)))) ## CHANGE TO LOGeUCLIDEAN IF NEEDED
        ##matImportances = mean(kernel.(sqrt.(distanceSqr.(logEuclidean, datmats1, refmats1))))
        matImportances = sum(kernel.(sqrt.(distanceSqr.(Fisher, datmats1, refmats1)), sigmaMat))

        Importances = push!(traitImportances, matImportances) 
        ## Importances = matImportances
        ##Importances = traitImportances 
        parameters = sigmasAll[j]
        weight = Importances
        push!(res, (parameters, weight))
        ## push!(res, Importances)
        next!(p)
    end
    res, unrun
end

# Perform simulation
Threads.@threads for i in 1:5 # Use number of available threads
    tree = trees[i]
    this_data = data[i]
    println("Running simulation for tree $i on thread $(threadid())")
tst, unrun = sim(5000, tree, this_data) ## CHANGE NUMBER OF PARTICLES IF NEEDED
##tst, unrun = sim(20000, tree, this_data) #final optimization
@save "C:/Users/suqt0/OneDrive/Desktop/UQ/JMMenura/paper_scripts/importance_simulated_scripts/simulated/Test Data sets and trees/newSim/BMmodelled$i.BMtestdata.s2.5_5.r4.sP10.logE.mean.jld2" tst unrun
    println("Simulation for tree $i completed and saved.")
end
