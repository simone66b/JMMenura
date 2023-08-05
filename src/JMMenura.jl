module JMMenura

#include("Basic_Functions.jl")
include("Basic_Functions_exper.jl")
include("Basic_Functions_each.jl")
include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
include("Tree_Modifying_Functions.jl")
include("Plotting.jl")
include("JMM_Bayesian.jl")
export menura_sim, menura_sim_exper, menura_sim_exper_mat_OU, menura_sim_exper_mat_wiener_process

export menura_sim_mat_OU_each

export trait_drift, trait_drift_brownian_motion, trait_diff, trait_diff_beta, trait_diff_cox_ingersoll_ross_gamma

export covariance_mat_drift, covariance_mat_diffusion

export apply_trait, apply_prior,apply_prior_descend, apply_trait_descend, get_trait

export plot_labelled

export Menura_bayesian

end # module 
