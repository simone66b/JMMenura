using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

n = 3

a1=1.0;
b1= 1.0;
tree1 = Ultrametric(6);
tree = rand(tree1); 
time_tot = 1.0;
tspan = (0.0, time_tot);

P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) 

alpha1 = repeat([1.0], n);
mu1 = repeat([0.0], n); ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n);
parms= (alpha=alpha1, sigma=sigma1)

reference_data = menura_sim_exper(alpha1, sigma1, mu1, P0, a1, b1, tree)
known_var = nothing

prior = [repeat([Normal(1,1)], n)..., repeat([Exponential(1)], n)...,
            repeat([Normal(0,1)], n)..., Wishart(100, Matrix(1I, n, n))]
apply_prior_descend(tree, prior, 1)

Menura_bayesian(reference_data[2], tree, known_var; max_iter = convert(Int, 2e4), x0 = nothing, 
    threshold = 30.0, n_particles= 200)



# for tip in getleaves(reference_data[1])
#     tip.data["trace"][end]
# end

# [tip.data["trace"][end] for tip in getleaves(reference_data[1])]
