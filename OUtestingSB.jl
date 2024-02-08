using DifferentialEquations, Distributions, Random, Plots, LinearAlgebra, PosDefManifold
pyplot()

G0 =  [1.0 0.8;
0.8 1.0]
##**1.0I(2)

muG = [3.00 0.8;
       0.8 2.0]
alpha = 0.0
sigma = 1.0I(2)
deltat = 0.001
its = 1000
sigmamat = fill(1/sqrt(2), 2, 2)

temp = []
push!(temp, G0)

W = WienerProcess(0.0, 0.0)
prob=NoiseProblem(W, (0.0, 1.0))
Wt1 = solve(prob; dt=0.001).u[1:its]
Wt2 = solve(prob; dt=0.001).u[1:its]
Wt3 = solve(prob; dt=0.001).u[1:its]
Wts = []
push!(Wts, zeros(2,2))
for i = 2:its
    Wt = [Wt1[i] Wt3[i]/sqrt(2);
          Wt3[i]/sqrt(2) Wt2[i]]  
    G = temp[i-1]
    sqrtG = sqrt(G)
    invsqrtG = inv(sqrtG)
    ans = sqrtG * exp(alpha * invsqrtG * muG * invsqrtG * deltat + 
    sigma * (sqrt(deltat) * (Wt-Wts[i-1]))) * sqrtG
    push!(Wts, Wt)
    push!(temp, ans)
end # i iterator

## temp is now a vector of 3 x 3 Hermitian matrices, sampled using the Mai's OU simulation.
vals11 = map(x -> x[1,1], temp)
vals12 = map(x -> x[1,2], temp) 
vals21 = map(x -> x[2,1], temp) ## should be same as vals12
vals22 = map(x -> x[2,2], temp)

plot(vals11)
plot!(vals12)
plot!(vals22)
plot!(vals21)
