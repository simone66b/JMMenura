module JMMenura

using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, KissABC
using PosDefManifold, GpABC, StatsPlots, ProgressBars

include("JMMABCparameters.jl")
include("Basic_Functions_exper.jl")
include("Basic_Functions_each.jl")
include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
include("Tree_Modifying_Functions.jl")
include("Plotting.jl")
include("JMM_Bayesian.jl")
include("Sim_Functions.jl")


export menura_sim, menura_sim_exper_mat_isospectral, menura_sim_exper_mat_OU, menura_sim_exper_mat_wiener_process

export menura_sim_mat_OU_each

export trait_drift, trait_drift_brownian_motion, trait_diff, trait_diff_beta, trait_diff_cox_ingersoll_ross_gamma

export covariance_mat_drift, covariance_mat_diffusion

export apply_trait, apply_prior,apply_prior_descend, apply_trait_descend, get_trait

export plot_labelled, plot_data, plot_g_mat_evol, plot_traits_cov, animate_data

export menura_bayesian, trait_mat_distance, get_data, menura_bayesian_trait_alpha, menura_bayesian_mat_alpha, menura_bayesian_trait_mat_alpha

export test_threshold

export mat_evol, trait_evol, mat_evol_skew_symmetric

export menura_para_descend!, menura_redone!

export get_priors, assemble_mat_parameters, assemble_trait_parameters, JMMABCAllSingleConstnat, JMMABCAlphaEqualConstant, JMMABCAlphaConstantEqual

end # module 
