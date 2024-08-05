######################
# Struct definitions #
######################

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
    mat_sigma_known::Number
    size::Int
end

"""
Stores the priors for a JMM ABC simulation where only alpha parameters and sigma parameter for matrix evolution are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is constant for all variables
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCAlphaEqualConstantSigma <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_prior::ContinuousUnivariateDistribution
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
    mat_sigma_known::Number
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
struct JMMABCAlphaSigmaDifferentConstant <: JMMABCparameters
    trait_alpha_prior::Array{<:ContinuousUnivariateDistribution}
    trait_mu_known::Array{<:Number}
    trait_sigma_prior::Array{<:ContinuousUnivariateDistribution}
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
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
struct JMMABCAlphaEqualConstantTraitBrownian <: JMMABCparameters
    trait_mu_known::Array{<:Number}
    trait_sigma_prior::ContinuousUnivariateDistribution
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_known::Number
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
Stores the priors for a JMM ABC simulation where only alpha parameters are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and can be different for different branches
- Other parameters can be different on the branches

NOTE: This setup requires knowing the internal names of the link nodes in a tree. These can be used using
the plot_labelled function included in this package.
Every dictionary must contain all modified nodes

Inputs:
TODO
"""
struct JMMABCAlphaDifferentEqual <: JMMABCparameters
    trait_alpha_prior::Dict{Int64, Any}
    trait_mu_known::Dict{Int64, Any}
    trait_sigma_known::Dict{Int64, Any}
    mat_alpha_prior::Dict{Int64, Any}
    mat_mu_known::Dict{Int64, Any}
    mat_sigma_known::Dict{Int64, Any}
    size::Int
end



"""
Stores the priors for a JMM ABC simulation using the isospectral model where only the alpha parameter is unknown.

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

"""
Stores the priors for a JMM ABC simulation using the isospectral model where only the alpha parameter is unknown.

Conditions:
- Alpha is unknown traits and is constant for all variables and branches
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCIsospectralAlphaAB <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_a_prior::ContinuousUnivariateDistribution
    mat_b_prior::ContinuousUnivariateDistribution
    size::Int
end

"""
Stores the priors for a JMM ABC simulation using the isospectral model where only the alpha parameter is unknown.

Conditions:
- Alpha is unknown traits and is constant for all variables and branches
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCIsospectralAlphaABTraitBrownian <: JMMABCparameters
    trait_mu_known::Array{<:Number}
    trait_sigma_prior::ContinuousUnivariateDistribution
    mat_a_prior::ContinuousUnivariateDistribution
    mat_b_prior::ContinuousUnivariateDistribution
    size::Int
end


"""
Stores the priors for a JMM ABC simulation using the isospectral model where only the alpha parameter is unknown.

Conditions:
- Alpha is unknown traits and is constant for all variables and branches
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCIsospectralAlphaABTraitOUDiff <: JMMABCparameters
    trait_alpha_prior::Array{<:ContinuousUnivariateDistribution}
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_a_prior::ContinuousUnivariateDistribution
    mat_b_prior::ContinuousUnivariateDistribution
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
struct JMMABCBrownian <: JMMABCparameters
    trait_alpha_prior::ContinuousUnivariateDistribution
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_mu_known::Array{<:Number}
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
struct JMMABCBrownianTraitsBrownian <: JMMABCparameters
    trait_mu_known::Array{<:Number}
    trait_sigma_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
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
struct JMMABCBrownianTraitsOUDiff <: JMMABCparameters
    trait_alpha_prior::Array{<:ContinuousUnivariateDistribution}
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    mat_mu_known::Array{<:Number}
    mat_sigma_prior::ContinuousUnivariateDistribution
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
struct JMMABCAlphaDifferentConstantMatrix <: JMMABCparameters
    trait_alpha_prior::Array{<:ContinuousUnivariateDistribution}
    trait_mu_known::Array{<:Number}
    trait_sigma_known::Array{<:Number}
    size::Int
end


##################
# No Trait types #
##################

"""
Stores the priors for a JMM ABC simulation without evolving traits where only alpha parameters are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is constant for all variables
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCAlphaEqualConstantNoTrait <: JMMABCparameters
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_known::Number
    size::Int
end

"""
Stores the priors for a JMM ABC simulation without evolving traits where only alpha and sigma parameters are unknown.

Conditions:
- Alpha is unknown for matrixes
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCAlphaEqualConstantNoTraitSigma <: JMMABCparameters
    mat_alpha_prior::ContinuousUnivariateDistribution
    mat_mu_known::Array{<:Number}
    mat_sigma_prior::ContinuousUnivariateDistribution
    size::Int
end

"""
Stores the priors for a JMM ABC simulation using the isospectral model without evolving traits where only the alpha parameter is unknown.

Conditions:
- Alpha is unknown traits and is constant for all variables and branches
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCIsospectralAlphaABNoTrait <: JMMABCparameters
    mat_a_prior::ContinuousUnivariateDistribution
    mat_b_prior::ContinuousUnivariateDistribution
    size::Int
end

"""
Stores the priors for a JMM ABC simulation without traits where only alpha parameter are unknown.

Conditions:
- Alpha is unknown for matrixes and traits and is constant for all variables
- Parameters are constant on all branches of a tree

Inputs:
TODO
"""
struct JMMABCBrownianNoTrait <: JMMABCparameters
    mat_mu_known::Array{<:Number}
    mat_sigma_prior::ContinuousUnivariateDistribution
    size::Int
end

############################################
# Functions to extract prior distributions #
############################################

"""
Returns the priors as a vector of continuous univariate distributions

Inputs: 
- parameters - An object which is a subtype of the JMMABCparameters type 

Outputs: 

A vector with elements of type ContinuousUnivariateDistribution
"""
function get_priors(parameters::JMMABCparameters) end

function get_priors(parameters::JMMABCAllEqualConstant)
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.trait_mu_prior, parameters.trait_sigma_prior, 
                parameters.mat_alpha_prior, parameters.mat_mu_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCAlphaEqualConstant) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaEqualConstantSigma) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_alpha_prior, 
                                                        parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCAlphaDifferentConstant) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior..., parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaSigmaDifferentConstant) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior..., parameters.trait_sigma_prior..., parameters.mat_alpha_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCAlphaEqualConstantTraitBrownian) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_sigma_prior, parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaConstantEqual) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaDifferentEqual) 
    trait_indexes = sort(collect(keys(parameters.trait_alpha_prior)), rev = true)
    mat_indexes = sort(collect(keys(parameters.mat_alpha_prior)), rev = true)
    return  Vector{ContinuousUnivariateDistribution}([[parameters.trait_alpha_prior[i] for i in trait_indexes]..., 
    [parameters.mat_alpha_prior[i] for i in mat_indexes]...])
end

function get_priors(parameters::JMMABCIsospectralAlpha) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior])
end

function get_priors(parameters::JMMABCIsospectralAlphaAB) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_a_prior, parameters.mat_b_prior])
end

function get_priors(parameters::JMMABCIsospectralAlphaABTraitBrownian) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_sigma_prior, parameters.mat_a_prior, parameters.mat_b_prior])
end

function get_priors(parameters::JMMABCIsospectralAlphaABTraitOUDiff) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior..., parameters.mat_a_prior, parameters.mat_b_prior])
end

function get_priors(parameters::JMMABCBrownian) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCBrownianTraitsBrownian) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_sigma_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCBrownianTraitsOUDiff) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.trait_alpha_prior..., parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCAlphaDifferentConstantMatrix) 
    return  Vector{ContinuousUnivariateDistribution}(parameters.trait_alpha_prior)
end

######################
# No Trait Functions #
######################

function get_priors(parameters::JMMABCAlphaEqualConstantNoTrait) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.mat_alpha_prior])
end

function get_priors(parameters::JMMABCAlphaEqualConstantNoTraitSigma) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.mat_alpha_prior, parameters.mat_sigma_prior])
end

function get_priors(parameters::JMMABCIsospectralAlphaABNoTrait) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.mat_a_prior, parameters.mat_b_prior])
end

function get_priors(parameters::JMMABCBrownianNoTrait) 
    return  Vector{ContinuousUnivariateDistribution}([parameters.mat_sigma_prior])
end



#####################################################################################
# Functions to assemble the parameters for a JMMenura simulation  using an OU model #
#####################################################################################

"""
Combines a random result from the relevant priors for a distribution with the constant parameters for evolving the traits

Inputs:
- parameters - An object which is a subtype of the JMMABCparameters type 
- prior_results - A vector containing random results of the priors specified by parameters

Outputs:

A named tuple with elements being the relevant parameters for the model specified by parameters
"""
function assemble_trait_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_trait_parameters(parameters::JMMABCAllEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = prior_results[2]*ones(parameters.size), sigma = prior_results[3]*ones(parameters.size)) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstantSigma, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(parameters.size)], mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaSigmaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(parameters.size)], mu = parameters.trait_mu_known, sigma = prior_results[(1 + parameters.size):(2 * parameters.size)]) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstantTraitBrownian, prior_results::Vector{<:Number}) 
    return (alpha = zeros(parameters.size), mu = parameters.trait_mu_known, sigma = prior_results[1]*ones(parameters.size)) 
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

function assemble_trait_parameters(parameters::JMMABCAlphaDifferentEqual, prior_results::Vector{<:Number}) 
    combined = Dict()
    for (i, ind) in enumerate(sort(collect(keys(parameters.trait_mu_known)), rev = true))
        combined[ind] = (alpha = prior_results[i]*ones(parameters.size)
        , mu = parameters.trait_mu_known[ind]
        , sigma = parameters.trait_sigma_known[ind])
    end
    return combined
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlpha, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlphaAB, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlphaABTraitBrownian, prior_results::Vector{<:Number}) 
    return (alpha = zeros(parameters.size), mu = parameters.trait_mu_known, sigma = prior_results[1]*ones(parameters.size)) 
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlphaABTraitOUDiff, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(parameters.size)], mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCBrownian, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1]*ones(parameters.size), mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCBrownianTraitsBrownian, prior_results::Vector{<:Number}) 
    return (alpha = zeros(parameters.size), mu = parameters.trait_mu_known, sigma = prior_results[1]*ones(parameters.size)) 
end

function assemble_trait_parameters(parameters::JMMABCBrownianTraitsOUDiff, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(parameters.size)], mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end

function assemble_trait_parameters(parameters::JMMABCAlphaDifferentConstantMatrix, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1:(parameters.size)], mu = parameters.trait_mu_known, sigma = parameters.trait_sigma_known) 
end


############
# No Trait #
############

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstantNoTrait, prior_results::Vector{<:Number}) 
    return nothing 
end

function assemble_trait_parameters(parameters::JMMABCAlphaEqualConstantNoTraitSigma, prior_results::Vector{<:Number}) 
    return nothing 
end

function assemble_trait_parameters(parameters::JMMABCIsospectralAlphaABNoTrait, prior_results::Vector{<:Number}) 
    return nothing
end

function assemble_trait_parameters(parameters::JMMABCBrownianNoTrait, prior_results::Vector{<:Number}) 
    return nothing 
end

"""
Combines a random result from the relevant priors for a distribution with the constant parameters for evolving the G matrix

Inputs:
- parameters - An object which is a subtype of the JMMABCparameters type 
- prior_results - A vector containing random results of the priors specified by parameters

Outputs:

A named tuple with elements being the relevant parameters for the model specified by parameters
"""
function assemble_mat_parameters(parameters::JMMABCparameters, prior_results::Vector{<:Number}) end

function assemble_mat_parameters(parameters::JMMABCAllEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[4], mu = prior_results[5]*ones(parameters.size, parameters.size), sigma = prior_results[6]) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[2], mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstantSigma, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[2], mu = parameters.mat_mu_known, sigma = prior_results[3]) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[end], mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaSigmaDifferentConstant, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[end-1], mu = parameters.mat_mu_known, sigma = prior_results[end]) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstantTraitBrownian, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[2], mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
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

function assemble_mat_parameters(parameters::JMMABCAlphaDifferentEqual, prior_results::Vector{<:Number})
    combined = Dict()
    offset = length(collect(keys(parameters.trait_mu_known)))
    for (i, ind) in enumerate(sort(collect(keys(parameters.mat_mu_known)), rev = true))
        combined[ind] = (alpha = prior_results[offset + i]
        , mu = parameters.mat_mu_known[ind]
        , sigma = parameters.mat_sigma_known[ind])
    end
    return combined
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlpha, prior_results::Vector{<:Number}) 
    return (a = parameters.mat_a,b = parameters.mat_b) 
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlphaAB, prior_results::Vector{<:Number}) 
    return (a = prior_results[2],b = prior_results[3]) 
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlphaABTraitBrownian, prior_results::Vector{<:Number}) 
    return (a = prior_results[2],b = prior_results[3]) 
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlphaABTraitOUDiff, prior_results::Vector{<:Number}) 
    return (a = prior_results[parameters.size + 1],b = prior_results[parameters.size + 2]) 
end


function assemble_mat_parameters(parameters::JMMABCBrownian, prior_results::Vector{<:Number}) 
    return (alpha = 0, mu = parameters.mat_mu_known, sigma = prior_results[2]) 
end

function assemble_mat_parameters(parameters::JMMABCBrownianTraitsBrownian, prior_results::Vector{<:Number}) 
    return (alpha = 0, mu = parameters.mat_mu_known, sigma = prior_results[2]) 
end

function assemble_mat_parameters(parameters::JMMABCBrownianTraitsOUDiff, prior_results::Vector{<:Number}) 
    return (alpha = 0, mu = parameters.mat_mu_known, sigma = prior_results[parameters.size+1]) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaDifferentConstantMatrix, prior_results::Vector{<:Number}) 
    return (alpha = 1, sigma = 0) 
end

############
# No Trait #
############

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstantNoTrait, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1], mu = parameters.mat_mu_known, sigma = parameters.mat_sigma_known) 
end

function assemble_mat_parameters(parameters::JMMABCAlphaEqualConstantNoTraitSigma, prior_results::Vector{<:Number}) 
    return (alpha = prior_results[1], mu = parameters.mat_mu_known, sigma = prior_results[2]) 
end

function assemble_mat_parameters(parameters::JMMABCIsospectralAlphaABNoTrait, prior_results::Vector{<:Number}) 
    return (a = prior_results[1],b = prior_results[2]) 
end

function assemble_mat_parameters(parameters::JMMABCBrownianNoTrait, prior_results::Vector{<:Number}) 
    return (alpha = zeros(parameters.size, parameters.size), mu = parameters.mat_mu_known, sigma = prior_results[1]) 
end

#########################################
# Bayesian simulation function creators #
#########################################

"""
Creates a function which can be used to perform approximate bayesian computation in the style given by GpABC

Inputs:
- tree - An object of type Tree which the simulation is to be performed on
- JMMpara - An object which is a subtype of the JMMABCparameters type 
- trait0 - A vector of Numbers which represents the intial values for the traits
- mat0 - A matrix of Numbers which represents the intial values for the G matrix

Optional inputs: 
- t0 - The inital time to start the simulation from
- each - If false the G matrix matrix is kept constant when evolving the traits along each of the branches. If true the G matrix is updated at each time step when evolving the traits
- dt - The time step to be used

Outputs:

A function which can be used in the framework provided by GpABC
"""
function create_bayesian_sim(tree, JMMpara::JMMABCparameters, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001) end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstant, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstantSigma, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaDifferentConstant, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()

        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaSigmaDifferentConstant, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()

        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstantTraitBrownian, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaConstantEqual, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = assemble_trait_parameters(JMMpara, parameter)

        mat_para = assemble_mat_parameters(JMMpara, parameter)

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaDifferentEqual, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = assemble_trait_parameters(JMMpara, parameter)

        mat_para = assemble_mat_parameters(JMMpara, parameter)

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end


function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlpha, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_isospectral(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end


function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlphaAB, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_isospectral(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlphaABTraitBrownian, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_isospectral(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlphaABTraitOUDiff, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_isospectral(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCBrownian, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end


function create_bayesian_sim(tree, JMMpara::JMMABCBrownianTraitsBrownian, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end


function create_bayesian_sim(tree, JMMpara::JMMABCBrownianTraitsOUDiff, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaDifferentConstantMatrix, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => assemble_trait_parameters(JMMpara, parameter)) 

        mat_para = Dict(tree.nodedict[root.name] => (assemble_mat_parameters(JMMpara, parameter)..., mu = mat0))

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(dt = dt), mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()

        return get_data_no_mat(sim)
    end
end


#############
# No Traits #
#############


function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstantNoTrait, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => nothing) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, nothing, mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data_no_trait(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCAlphaEqualConstantNoTraitSigma, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => nothing) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, nothing, mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data_no_trait(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCIsospectralAlphaABNoTrait, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => nothing) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, nothing, mat_evol_isospectral(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data_no_trait(sim)
    end
end

function create_bayesian_sim(tree, JMMpara::JMMABCBrownianNoTrait, trait0, mat0; t0 = 0.0, each = true, 
    dt = 0.001)
    function bayesian_menura!(parameter)

        root = getroot(tree)

        trait_para = Dict(tree.nodedict[root.name] => nothing) 

        mat_para = Dict(tree.nodedict[root.name] => assemble_mat_parameters(JMMpara, parameter)) 

        sim = menura_parameter_descend!(mat_para, trait_para, tree, nothing, mat_evol_affine(dt = dt), t0, trait0, mat0, each)
        
        GC.gc()
        
        return get_data_no_trait(sim)
    end
end