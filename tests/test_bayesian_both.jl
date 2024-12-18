using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances
pyplot()


Random.seed!(1)

# number of parameters
n = 3

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))

# traits needed to evolve traits
alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n)
sigma1 = repeat([1.0], n)

# create trait dictionary
trait_parameters_true = Dict(11 => (alpha = alpha1, mu = mu1, sigma = sigma1))
trait_parameters = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha = 1 .* Matrix(1I, n, n)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters_true = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))
mat_parameters = (mu = mat_mu, sigma = mat_sigma)

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

α_priors = [Gamma(2,0.25), Gamma(2,0.25)]

ref_sim = menura_para_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

ref_data = get_data(ref_sim[1])


@time out = menura_bayesian_trait_mat_alpha(ref_data, tree1, trait_parameters
    ,mat_parameters, α_priors, mu1, P0, 300, n, 300.0)

scatter([out.population[i][1] for i in 1:length(out.population)], out.distances)

heatmap(out.population, out.distances)
