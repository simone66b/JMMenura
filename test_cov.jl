using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, DifferentialEquations

pyplot()

Random.seed!(5)

n = 6

a1=1.0;
b1= 1.0;
tree1 = Ultrametric(20);
    tree = rand(tree1); 
       time_tot = 1.0;
    tspan = (0.0, time_tot);

P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) # Change to Wishart Distribution

alpha1 = repeat([1.0], n);
mu1 = repeat([0.0], n); ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n);
parms= (alpha=alpha1, sigma=sigma1)

u0=zeros(size(P0))

p1 = (alpha=alpha1, sigma=sigma1, mu=mu1, mat=P0, a=a1, b=b1)
dt = 0.001


lowertri = LowerTriangular(P0)
uppertri = - UpperTriangular(P0)
skewsymm = lowertri + uppertri
W = WienerProcess(0.0, 0.0, 0.0)

pp = (A=skewsymm, B=skewsymm, a=p1.a, b=p1.b) ## skew symmetric matrices not necessarily the same.
prob = SDEProblem(covariance_mat_drift, covariance_mat_diffusion, u0, tspan, p=pp, noise=W,
                noise_rate_prototype=zeros(size(P0))) ## setup SDE problem
msol = solve(prob, EM(), p=pp, dt=dt)
Omegas = exp.(msol.u)

Ps = [Omegai * P0 * Omegai' for Omegai in Omegas]

upper = maximum(maximum.(Ps))
lower = minimum(minimum.(Ps))

anim = @animate for Pi in Ps
    heatmap(Pi, yflip = true, clims = (lower, upper))
end every 10

gif(anim, "cov_evolution_test.gif", fps=30)

nodes = [getancestors(exampledat1[1], exampledat1[1].nodes[3])..., exampledat1[1].nodes[3]]

anim = @animate for (i, node) in enumerate(nodes)
    heatmap(cor(convert(Matrix{Float64},node.data["matrix"])), title = "$i", yflip = true, clims = (-1, 1))
end 

gif(anim, "Cov_evol.gif", fps = 0.5)

cov_mats = [convert(Matrix{Float64},node.data["matrix"]) for node in nodes]

lower = minimum([minimum(cov_mat) for cov_mat in cov_mats])
upper = maximum([maximum(cov_mat) for cov_mat in cov_mats])

anim = @animate for (i, cov_mat) in enumerate(cov_mats)
    heatmap(cov_mat, title = "$i", yflip = true, clims = (lower, upper))
end 

gif(anim, "Cov_evol.gif", fps = 0.5)