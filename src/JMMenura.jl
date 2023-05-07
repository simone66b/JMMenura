module JMMenura

#include("Basic_Functions.jl")
include("Basic_Functions_exper.jl")
include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
export menura_sim, menura_sim_exper

export trait_drift, trait_drift_brownian_motion, trait_diff, trait_diff_beta, trait_diff_cox_ingersoll_ross_gamma

end # module 
