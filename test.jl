using .JMMenura
using Phylo

a1=1.0;
b1= 1.0;
tree1 = Ultrametric(20);
    tree = rand(tree1); 
       time_tot = 1.0;
    tspan = (0.0, time_tot);

P0 = [0.329 0.094 -0.083 -0.089 0.293 0.079 0.208 0.268;
0.094 0.449 0.349 0.24 0.071 0.075 0.03 0.009;
-0.083 0.349 1.426 0.487 -0.371 -0.098 -0.053 -0.172;
-0.089 0.24 0.487 0.546 -0.168 0.017 -0.051 -0.081;
0.293 0.071 -0.371 -0.168 1.441 1.008 0.904 0.945;
0.079 0.075 -0.098 0.017 1.008 1.087 0.731 0.78;
0.208 0.03 -0.053 -0.051 0.904 0.731 0.809 0.783;
0.268 0.009 -0.172 -0.081 0.945 0.78 0.783 0.949];

alpha1 = repeat([1.0], 8);
mu1 = repeat([0.0], 8); ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], 8);
parms= (alpha=alpha1, sigma=sigma1)

menura_sim(alpha1, sigma1, mu1, P0, a1, b1, tree)
