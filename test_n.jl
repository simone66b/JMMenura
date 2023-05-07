using .JMMenura
using Phylo, Distributions, Random

Random.seed!(1)

n = 9

a1=1.0;
b1= 1.0;
tree1 = Ultrametric(20);
    tree = rand(tree1); 
       time_tot = 1.0;
    tspan = (0.0, time_tot);

P0 = rand(Uniform(-1,1), n,n)

alpha1 = repeat([1.0], n);
mu1 = repeat([0.0], n); ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n);
parms= (alpha=alpha1, sigma=sigma1)


@time menura_sim(alpha1, sigma1, mu1, P0, a1, b1, tree)


@time menura_sim_exper(alpha1, sigma1, mu1, P0, a1, b1, tree, x0 = [1,1,1,1,1,1,1,1,1])
