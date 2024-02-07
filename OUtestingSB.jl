using Distributions, Random, Plots, LinearAlgebra
pyplot()

G0 =  [0.2919173 0.2128847 0.0000000;
       0.2128847 1.0333177 0.0000000;
       0.0000000 0.0000000 0.8134763]

muG = [1.4935775 0.5337943 0.000000;
       0.5337943 0.7247179 0.000000;
       0.0000000 0.0000000 0.630625]

sigmamat = fill(1/sqrt(2), size(G0))
sigmamat[diagind(sigmamat)] .= 1.0
alpha = 5
sigma = .01 .*I(3)
deltat = 0.001
its = 1000

temp = []
push!(temp, G0)

for i = 2:its
    Wt =  randn(9)
    Wt = sigmamat .* Hermitian(reshape(Wt, (3,3)))
    G = temp[i-1]
    sqrtG = sqrt(G)
    invsqrtG = inv(sqrtG)
    ans = sqrtG * exp(alpha * invsqrtG * muG * invsqrtG* deltat + sigma * sqrt(deltat)* Wt) * sqrtG
    push!(temp, ans)
end # i iterator
plot(map(x -> x[1,1], temp))
## temp is now a vector of 3 x 3 Hermitian matrices, sampled using the Mai's OU simulation.
