using GpABC
using Distributions
using DifferentialEquations
using Plots

sim_result = SimulatedABCRejection(reference_data, simulator_function, priors, threshold, n_particles;
    max_iter=convert(Int, 2e6),
    write_progress=true)

"""
A function which returns the required priors from a tree

Could make it so you can choose which ones
"""
function JMM_get_priors(tree)
    nothing 
end



function Menura_bayesian(reference_data, tree, known_var; trait_drift = trait_drift, trait_diff = trait_diff,
    matrix_drift = covariance_mat_drift, max_iter = convert(Int, 2e6))
    #For unknown functions get variables
    priors = JMM_get_priors(tree)

    # Define function in here - more customisable
    # Defines a version of the simulation function that only takes one input 
    function bayesian_menura_sim(paramaters)
    end
    SimulatedABCRejection(reference_data, bayesian_menura_sim, priors, threhold, n_particles, max_iter = max_iter, write_progress = true)
end
