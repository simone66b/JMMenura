# Creates a basic evolution plot
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, StatsPlots, MultivariateStats, Plots.PlotMeasures

pyplot()

Random.seed!(1)

# number of parameters
n = 7

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
trait_parameters = Dict(11 => (alpha = alpha1, mu = mu1, sigma = sigma1))

# Variables needed for OU matrix model
mat_alpha = 1 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

menura_para_descend!(mat_parameters, trait_parameters, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)


p1 = plot(tree1)

leaves = getleaves(tree1)


l = @layout [ a{0.7w} [grid(6, 1)]]
p1 = plot([plot(tree1), [plot() for _ in 1:6]...]..., layout = l)

for (i, leaf) in enumerate(reverse(leaves))
    cov = leaf.data["mat_trace"][end]
    trait = leaf.data["trait_trace"][end]

    pca = fit(PCA, cov, maxoutdim = 2)

    pca_mu = predict(pca, trait)

    proj = projection(pca)

    # append!(covps, plot(covellipse(pca_mu,proj'*cov*proj)))
    p1 = covellipse!(pca_mu,proj'*cov*proj, subplot = i+1, legend = false, title = leaf.name)
end


display(p1)