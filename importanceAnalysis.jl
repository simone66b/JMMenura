using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions

cd("/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/Test Data sets and trees/")

files = readdir(".")
## jld2_files = filter(f -> endswith(f, "a.jld2"), files)

i = 1
## for i in 1:5
    pattern = ".*modelled$i\\.BMtestdata"
    re = Regex(pattern)
    matches = filter(f -> occursin(re, f), files)
    BMfile  = findfirst(f -> occursin("BM", f), matches)
    ISOfile = findfirst(f -> occursin("ISO", f), matches)
    OUfile  = findfirst(f -> occursin("OU", f), matches)

    BM = JLD2.load(matches[BMfile])
    ISO = JLD2.load(matches[ISOfile])
    OU = JLD2.load(matches[OUfile])

    BM = vcat(BM["tst"]...)
    BM = vcat(BM...)
    
    ISO = vcat(ISO["tst"]...)
    ISO = vcat(ISO...)

    OU = vcat(OU["tst"]...)
    OU = vcat(OU...)
## end

# lw = log.(weights)
# maxlw = maximum(lw)
# logZ = maxlw + log(sum(exp.(lw .- maxlw))) - log(N)