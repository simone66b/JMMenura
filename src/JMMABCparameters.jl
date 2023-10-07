
"""
Abstract Super type which for types which house the conditions and assumptions of various ABC models for JMM
"""
abstract type JMMABCparameters end 

"""
Stores the priors for a JMM ABC simulation where all parameters are unknown. Also assumes that each parameter is constant for all variables.

Inputs:
TODO
"""
struct JMMABCAllEqual <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_prior::ContinuousUnivariateDistribution
    trait_sigma_prior::ContinuousUnivariateDistribution
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_prior::ContinuousUnivariateDistribution
    mat_sigma_prior::ContinuousUnivariateDistribution
    size::Int
end

"""
Stores the priors for a JMM ABC simulation where only alpha parameters are unknown. Also assumes that each parameter is constant for all variables.

Inputs:
TODO
"""
struct JMMABCAlphaEqual <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Number
    trait_sigma_known::Number
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Number
    mat_sigma_known::Number
    size::Int
end

############################################
# Functions to extract prior distributions #
############################################

"""
Returns the priors as a vector of continuous univariate distributions
"""
function get_priors(parameters::JMMABCparameters) end

function get_priors(parameters::JMMABCAllEqual)
    return  [parameters.trait_alpha_prior, parameters.trait_mu_prior, parameters.trait_sigma_prior, 
                parameters.mat_alpha_prior, parameters.mat_mu_prior, parameters.mat_sigma_prior]
end

function get_priors(parameters::JMMABCAlphaEqual) 
    return  [parameters.trait_alpha_prior, parameters.mat_alpha_prior]
end

#####################################################################################
# Functions to assemble the parameters for a JMMenura simulation  using an OU model #
#####################################################################################

"""
"""
function assemble_trait_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_trait_parameters(parameters::JMMABCAllEqual, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = prior_results[2]*ones(parameters.size), sigma = prior_results[3]*ones(parameters.size)) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqual, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

"""
"""
function assemble_mat_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_mat_parameters(parameters::JMMABCAllEqual, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[4]*ones(parameters.size, parameters.size), mu = prior_results[5]*ones(parameters.size, parameters.size), sigma = prior_results[6]*ones(parameters.size, parameters.size)) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqual, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[2]*ones(parameters.size, parameters.size), mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end