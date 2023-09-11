using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, Distances
pyplot()

# Checking trait application

function menura_bayesian_trait_alpha_no_sim(reference_data, tree, trait_known_para, alpha_prior
    ,mat_known_para, trait0, mat0, max_iter, length)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)
    
    trait_para = Dict(11 => (trait_known_para..., alpha = parameter*ones(length))) 

    sim = menura_para_descend!(mat_known_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end
end

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
trait_parameterstrue = Dict(11 => (alpha = alpha1, mu = mu1, sigma = sigma1))
trait_parameters1 = (mu = mu1, sigma = sigma1)

# Variables needed for OU matrix model
mat_alpha = 1 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.1) .* ones(n,n)
mat_mu = copy(P0)

# create matrix dictionary
mat_parameters = Dict(11 => (alpha = mat_alpha, mu = mat_mu, sigma = mat_sigma))

mat_evol_func = mat_evol()
trait_evol_func = trait_evol()

α_prior = Uniform(0,1)

ref_sim = menura_para_descend!(mat_parameters, trait_parameterstrue, tree1, trait_evol_func, mat_evol_func, 0.0, mu1, P0, true)

ref_data = get_data(ref_sim[1])


@time out = menura_bayesian_trait_alpha(ref_data, tree1, trait_parameters1, α_prior
    ,mat_parameters, mu1, P0, 50000, n)

histogram(out.population, normalize = :pdf, legend = false)