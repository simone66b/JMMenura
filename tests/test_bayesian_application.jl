using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances


Random.seed!(1)

# number of parameters
n = 3

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# inital parameters
# G matrixes
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))
P1 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))

# trait parameters
t_alpha = Gamma(2, 0.25)
t_mu1 = repeat([0.0], n)
t_sigma1 = repeat([1.0], n)

t_mu2 = repeat([2.0], n)
t_sigma2 = repeat([2.0], n)

t_mu = Dict(11 => t_mu1, 5 => t_mu2)
t_sigma = Dict(11 => t_sigma1, 5 => t_sigma2)

# matrix parameters
m_alpha = Gamma(2, 0.25)
m_mu1 = copy(P0)
m_sigma1 = 1.0*ones(n,n)

m_mu2 = copy(P1)
m_sigma2 = 2.0*ones(n,n)

m_mu = Dict(11 => m_mu1, 5 => m_mu2)
m_sigma = Dict(11 => m_sigma1, 5 => m_sigma2)

# create parameter struct
a = JMMABCAlphaConstantEqual(t_alpha, t_mu, t_sigma, m_alpha, m_mu, m_sigma, n)

# set up evolution Functions
mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

# create reference data
mat_parameters_true = assemble_mat_parameters(a, [0.0, 1.0])

trait_parameters_true = assemble_trait_parameters(a, [0.0, 1.0])

ref_sim = menura_para_descend!(mat_parameters_true, trait_parameters_true, tree1, trait_evol_func, mat_evol_func, 0.0, t_mu1, P0, false)

ref_data = get_data(ref_sim[1])

# test for required thresholds

thresholds = test_threshold(ref_data, tree1, a, t_mu1, P0, 100)

histogram(thresholds)

@time out = menura_bayesian(ref_data, tree1, parameters, mu1, P0, 175.0, 20)

