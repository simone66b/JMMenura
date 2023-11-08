
"""
Abstract Super type which for types which house the conditions and assumptions of various ABC models for JMM
"""
abstract type JMMABCparameters end 

"""
Stores the priors for a JMM ABC simulation where all parameters are unknown. Also assumes that each parameter is constant for all variables.

Inputs:
TODO
"""
struct JMMABCAllEqualConstant <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_prior::ContinuousUnivariateDistribution
    trait_sigma_prior::ContinuousUnivariateDistribution
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_prior::ContinuousUnivariateDistribution
    mat_sigma_prior::ContinuousUnivariateDistribution
    size::Int
end

"""
Stores the priors for a JMM ABC simulation where only alpha parameters are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is constant for all variables
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCAlphaEqualConstant <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_known::Array{<:Number}
    size::Int
end

"""
Stores the priors for a JMM ABC simulation where only alpha parameters are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is different for all variables
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCAlphaDifferentConstant <: JMMABCparameters
    trait_alpha_prior::Array{<:ContinuousUnivariateDistribution}
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_known::Array{<:Number}
    size::Int
end


"""
Stores the priors for a JMM ABC simulation where only alpha parameters are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is constant for all variables and branches
- Other parameters can be different on the branches

NOTE: This setup requires knowing the internal names of the link nodes in a tree. These can be used using
the plot_labelled function included in this package.
Every dictionary must contain all modified nodes

Inputs:
TODO
"""
struct JMMABCAlphaConstantEqual <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Dict{Int64, Any}
    trait_sigma_known::Dict{Int64, Any}
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Dict{Int64, Any}
    mat_sigma_known::Dict{Int64, Any}
    size::Int
end


"""
Stores the priors for a JMM ABC simulation using the isospectral model where only  the alpha parameter is unknown.

Conditions:
- Alpha is unknown traits and is constant for all variables and branches
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCIsospectralAlpha <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_a::Number
    mat_b::Number
    size::Int
end


############################################
# Functions to extract prior distributions #
############################################

"""
Returns the priors as a vector of continuous univariate distributions
"""
function get_priors(parameters::JMMABCparameters) end

function get_priors(parameters::JMMABCAllEqualConstant)
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.trait_mu_prior, parameters.trait_sigma_prior, 
                parameters.mat_alpha_prior, parameters.mat_mu_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCAlphaEqualConstant) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaDifferentConstant) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior..., parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaConstantEqual) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCIsospectralAlpha) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior])
end

#####################################################################################
# Functions to assemble the parameters for a JMMenura simulation  using an OU model #
#####################################################################################

"""
"""
function assemble_trait_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_trait_parameters(parameters::JMMABCAllEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = prior_results[2]*ones(parameters.size), sigma = prior_results[3]*ones(parameters.size)) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(length(prior_results)-1)], mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaConstantEqual, prior_results::Vector{<:Number}) 
    combined = Dict()
    for i in collect(keys(parameters.trait_mu_known))
        combined[i] = (alpha = prior_results[1]*ones(parameters.size)
        , mu = parameters.trait_mu_known[i]
        , sigma = parameters.trait_sigma_known[i])
    end
    return combined
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlpha, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

"""
"""
function assemble_mat_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_mat_parameters(parameters::JMMABCAllEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[4]*ones(parameters.size, parameters.size), mu = prior_results[5]*ones(parameters.size, parameters.size), sigma = prior_results[6]*ones(parameters.size, parameters.size)) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[2]*ones(parameters.size, parameters.size), mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[end]*ones(parameters.size, parameters.size), mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaConstantEqual, prior_results::Vector{<:Number})
    combined = Dict()
    for i in collect(keys(parameters.mat_mu_known))
        combined[i] = (alpha = prior_results[2]*ones(parameters.size, parameters.size)
        , mu = parameters.mat_mu_known[i]
        , sigma = parameters.mat_sigma_known[i])
    end
    return combined
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlpha, prior_results::Vector{<:Number}) 
    return (a = parameters.mat_a,b = parameters.mat_b) 
end


#########################################
# Bayesian simulation function creators #
#########################################

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstant, trait0, mat0; t0 = 0.0, each = false, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_para_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol(dt = dt), t0, trait0, mat0, each)

        return get_data(sim[1])
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaDifferentConstant, trait0, mat0; t0 = 0.0, each = false, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_para_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol(dt = dt), t0, trait0, mat0, each)

        return get_data(sim[1])
    end
end


function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlpha, trait0, mat0; t0 = 0.0, each = false, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_para_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_skew_symmetric(dt = dt), t0, trait0, mat0, each)

        return get_data(sim[1])
    end
end


