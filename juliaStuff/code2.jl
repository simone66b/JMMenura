using Pkg
Pkg.activate(".") ## from inside juliaStuff
Pkg.develop(path="/home/simoneb/Desktop/JMMenura/")
using JMMenura

using StatsBase, Plots, JLD2, Distributions
pyplot() 
@load "s_para_result.jld2"

mat = s_para_result[1].population[1]
wts = s_para_result[1].weights[1]/sum(s_para_result[1].weights[1]) ## do we really need to use this? Or just conclude that wts are very small so don't contribute much to the posterior? Not sure.
d = Truncated(Normal(0, 10), 0, Inf)
histogram(mat[:,4], normalize=true)
plot!(d)

