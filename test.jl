using DifferentialEquations, Phylo, Plots, Distributions, Distances, KissABC; pyplot();


tst1 = Factored(Normal(0,1), Normal(3,5), Normal(4,2))

nt  = ntuple(x->Normal(0, 3), 5)

tst1 = [Truncated(Normal(0,1), 0, Inf), Truncated(Normal(3,5), 0, Inf), Truncated(Normal(4,2), 0, Inf)]

Factored(repeat(tst1, 3)...,)
