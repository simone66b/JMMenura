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
    ans = sqrtG*exp((alpha*log(invsqrtG*muG*invsqrtG)*deltat + sigma*sqrt(deltat)*Wt)*invsqrtG)*sqrtG

    push!(Wts, Wt)
    push!(temp, ans)
end # i iterator

## temp is now a vector of 3 x 3 Hermitian matrices, sampled using the Mai's OU simulation.
vals11 = map(x -> x[1,1], tst)
vals12 = map(x -> x[1,2], tst) 
vals21 = map(x -> x[2,1], tst) ## should be same as vals12
vals22 = map(x -> x[2,2], tst)

plot(vals11)
plot!(vals12)
plot!(vals22)
plot!(vals21)


#############################################################################################

using Distributions, LinearAlgebra, PosDefManifold

function myfunc()
    alpha = 1
    sigma = [0.1 0 ; 0 0.1]
    mu = [2 -1 ; -1 2]
    dt = 0.0001
    G_0 = [1 0 ; 0 2]
    N = 100000

    Gs = [Matrix(1.0I, 2,2) for _ in 1:(N+1)]

    Gs[1] = G_0

    Sigma = fill(1.0 / sqrt(2.0), size(G_0, 1), size(_G0, 1));
    Sigma[diagind(Sigma)] .= 1.0;

    for i in 2:(N+1) 
        W_t_seed = rand(Normal(0,1), 3)
        W_t = [W_t_seed[1] W_t_seed[2]/sqrt(2) ; W_t_seed[2]/sqrt(2) W_t_seed[3]]
        sqrt_G = sqrt(Gs[i-1])
        inv_sqrt_G = I/sqrt_G
        # Gs[i] = sqrt_G*exp(inv_sqrt_G*(alpha*sqrt_G*log(inv_sqrt_G*mu*inv_sqrt_G)*sqrt_G*dt + sigma*sqrt(dt)*W_t)*inv_sqrt_G)*sqrt_G
        Gs[i] = sqrt_G*exp((alpha*log(inv_sqrt_G*mu*inv_sqrt_G)*dt + sigma*sqrt(dt)*W_t)*inv_sqrt_G)*sqrt_G

    end

    return Gs
end