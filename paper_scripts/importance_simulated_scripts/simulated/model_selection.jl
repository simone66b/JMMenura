using Pkg
cd("/Users/coope/OneDrive/Documents/Uni/Phylogenetics_coding/importance_scripts")
## cd("/home/simoneb/Desktop/JMMenura")
## Pkg.develop(path="/home/simoneb/Desktop/JMMenura")
using Phylo, Distributions, Pkg, Plots, DataFrames, XLSX, StatsBase, JLD2, LinearAlgebra, DifferentialEquations
using PosDefManifold, ProgressMeter, StatsPlots, Random, Distributions

##################################
# Load reference simulation data #
##################################
@load "./simulated/BM_model_larger_trait_3/BM_model_larger_trait_3_importance_result.jld2" tst
BM_result = tst
@load "./simulated/OU_diff_larger_trait_3/OU_diff_larger_trait_3_importance_result.jld2" tst
OU_result = tst
@load "./simulated/ISO_model_larger_trait_3/ISO_model_larger_trait_3_importance_result.jld2" tst
ISO_result = tst

wtsISO = [ISO_result[i][2] for i in 1:5000]
wtsBM = [BM_result[i][2] for i in 1:5000]
wtsOU = [OU_result[i][2] for i in 1:5000]

sum(sum(wtsISO))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO))) # 0.21084651323037476
sum(sum(wtsBM))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO))) # 0.1353594413829308
sum(sum(wtsOU))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO))) # 0.6537940453866945 