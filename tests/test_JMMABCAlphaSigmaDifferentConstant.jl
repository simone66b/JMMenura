using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances, StatsBase
pyplot()


Random.seed!(1)

# number of parameters
function main()
    n = 4

    # Creating tree
    tree = Ultrametric(7)
    Random.seed!(1)
    tree1 = rand(tree)
    time_tot = 1.0
    tspan = (0.0, time_tot)

    # G matrix
    P0 =rand(Wishart(100, Matrix(1I, n, n)  ))
    alpha1 = 1.5*rand()*ones(n)
    mu1 = repeat([0.0], n)
    sigma1 = repeat([sqrt(2)], n)

    # create trait dictionary
    trait_parameters_true = Dict(13 => (alpha = alpha1, mu = mu1, sigma = sigma1))
    trait_parameters = (mu = mu1)

    # Variables needed for OU matrix model
    mat_alpha1 = 1
    mat_sigma = sqrt(2)
    mat_mu = copy(P0)

    # create matrix dictionary
    mat_parameters_true = Dict(13 => (alpha = mat_alpha1, mu = mat_mu, sigma = mat_sigma))
    mat_parameters = (mu = mat_mu)

    mat_evol_func = mat_evol_affine(dt = 0.005)
    trait_evol_func = trait_evol(dt =0.005)

    α_prior = Gamma(2,0.25)
    sigma_prior = Gamma(2,0.25)
    ther_ref_sim = menura_parameter_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, ones(n), P0, true)

    return ther_ref_sim
end



plot_data(tree1, 1, ylim = (-2, 10), zlim = (-2.0, 2.0), legend = false, reuse = false)


ther_ref_data = get_all_data(ther_ref_sim)

parameters = JMMABCAlphaSigmaDifferentConstant(repeat([α_prior], n) , mu1,repeat([sigma_prior], n), α_prior , mat_mu, sigma_prior, n)

thresholds = test_threshold(ther_ref_data, tree1, parameters, mu1, P0, 20, dt = 0.005
, distance_function = trait_mat_distance(parameters.size,nnodes(tree1)), summary_function = get_all_data)

threshold = quantile(thresholds, (0:2)./length(thresholds))[2]

out = menura_bayesian(ther_ref_data, tree1, parameters, ones(n), P0, threshold, 10, dt = 0.005
, distance_function = trait_mat_distance(parameters.size,nnodes(tree1)), summary_function = get_all_data, max_iter = 4000)

@time out2 = menura_bayesian(ther_ref_data, tree1, parameters, ones(n), P0, threshold, 50, dt = 0.005)
