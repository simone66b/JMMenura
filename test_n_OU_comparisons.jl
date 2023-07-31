using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

n = 4

a1=1.0
b1= 1.0
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
Random.seed!(1)
tree2 = rand(tree)
Random.seed!(1)
tree3 = rand(tree)
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

@time exampledat1 = menura_sim_exper_mat_OU(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma, mat_mu, tree1)
@time exampledat2 = menura_sim_mat_OU_each(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma,
                                             mat_mu, tree2, small_dt_scale = 100)
@time exampledat3 = menura_sim_mat_OU_each(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma,
                                             mat_mu, tree3, small_dt_scale = 10000)

testnodesdata = [getnodes(exampledat1[1]), getnodes(exampledat2[1]), getnodes(exampledat3[1])]

for testnodes in testnodesdata
    plot1 = plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
    for i in testnodes
        u1= i.data["trace"];
        uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
        myt = i.data["timebase"];
        plot!(myt, uu1[:, 1], uu1[:,3]);
    end; # for
    current()
    display(plot1)
end

for testnodes in testnodesdata 
    plot2 = plot(xlim=(-2.0,2.0), ylim= (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
    for i in testnodes
        u1= i.data["trace"];
        uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
        myt = i.data["timebase"];
        plot!(uu1[:,3], uu1[:, 1], uu1[:,2]);
    end;
    current()
    display(plot2)
end
# for node in exampledat[1].nodes
#     h = heatmap(convert(Matrix,node.data["matrixes"][end]), yflip = true, clims = (-1, 1))
#     current()
#     display(h)
# end