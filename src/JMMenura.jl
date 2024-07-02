module JMMenura

using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots
using PosDefManifold, GpABC, StatsPlots, ProgressBars, StochasticDiffEq

include("Preallocation.jl")
include("JMMABCparameters.jl")
include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
include("Tree_Modifying_Functions.jl")
include("Plotting.jl")
include("JMM_Bayesian.jl")
include("Sim_Functions.jl")


# export menura_sim, menura_sim_exper_mat_isospectral, menura_sim_exper_mat_OU, menura_sim_exper_mat_wiener_process

# export menura_sim_mat_OU_each

export trait_drift_mean_reversion, trait_drift_brownian_motion, trait_diffusion_brownian_motion, 
        trait_diffusion_cox_ingersoll_ross_gamma, trait_diff_beta

export matrix_drift_isospectral, matrix_diffusion_isospectral, matrix_drift_mean_reversion, matrix_diffusion_brownian_motion

export apply_trait, apply_prior,apply_prior_descend, apply_trait_descend, get_trait

export plot_labelled, plot_data, plot_g_mat_evol, plot_traits_cov, animate_data

export menura_bayesian, trait_mat_distance, mat_distance, trait_distance, get_data, get_data_no_trait, get_data_no_mat

export test_threshold, trait_mat_distance_scaled

export mat_evol, trait_evol, mat_evol_isospectral, no_trait_evol, mat_evol_affine

export menura_parameter_descend!, menura!, preallocate_tree!

export get_priors, assemble_mat_parameters, assemble_trait_parameters, create_bayesian_sim

export JMMABCAllSingleConstant, JMMABCAlphaEqualConstant, JMMABCAlphaConstantEqual, JMMABCIsospectralAlpha, JMMABCAlphaDifferentConstant, JMMABCIsospectralAlphaAB, JMMABCBrownian

export JMMABCAlphaEqualConstantSigma, JMMABCAlphaDifferentConstantMatrix

export JMMABCAlphaEqualConstantNoTrait, JMMABCIsospectralAlphaABNoTrait, JMMABCBrownianNoTrait, JMMABCAlphaEqualConstantNoTraitSigma

export JMMABCAlphaEqualConstantTraitBrownian, JMMABCBrownianTraitsBrownian, JMMABCBrownianTraitsOUDiff, JMMABCIsospectralAlphaABTraitBrownian, JMMABCIsospectralAlphaABTraitOUDiff

end # module 
