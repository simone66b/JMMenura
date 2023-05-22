module JMMenura

#include("Basic_Functions.jl")
include("Basic_Functions_exper.jl")
include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
include("Tree_Modifying_Functions.jl")
include("Plotting.jl")
export menura_sim, menura_sim_exper

export trait_drift, trait_drift_brownian_motion, trait_diff, trait_diff_beta, trait_diff_cox_ingersoll_ross_gamma

export apply_trait, apply_prior,apply_prior_descend, apply_trait_descend

export plot_labelled

end # module 
