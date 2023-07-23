using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

n = 4

a1=1.0
b1= 1.0
tree1 = Ultrametric(6)
tree = rand(tree1)
time_tot = 1.0
tspan = (0.0, time_tot)

P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) # Change to Wishart Distribution

alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n) ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n)
parms= (alpha=alpha1, sigma=sigma1)

mat_alpha = 12.0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)



@time exampledat = menura_sim_exper_mat_OU(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma, mat_mu, tree)

p1 = plot(tree, size = (400, 800),
    ## markersize = 20, 
      series_annotations = text.(1:nnodes(tree), 15, :center, :center,
                                 :white),  linewidth=5, showtips=false)
 current()
display(p1)

testnodes = getnodes(exampledat[1]);

p2 = plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
    myt = i.data["timebase"];
    plot!(myt, uu1[:, 1], uu1[:,3]);
end; # for
current()
display(p2)

p3 = plot(xlim=(0.0,1.0), ylim= (-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
    myt = i.data["timebase"];
    plot!(myt, uu1[:, 3]);
end; # for
current()
display(p3)

p4 = plot(ylim=(-2.0, 2.0), xlim= (-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
    plot!(uu1[:,1], uu1[:, 3]);
end;
current()
display(p4)

p5 = plot(xlim=(-2.0,2.0), ylim= (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
    myt = i.data["timebase"];
    plot!(uu1[:,3], uu1[:, 1], uu1[:,2]);
end;
current()
display(p5)