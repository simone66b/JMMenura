using GpABC
using Distributions
using DifferentialEquations
using Plots


"""
A function which returns the required priors from a tree

Could make it so you can choose which ones
"""
function JMM_get_priors(tree, known_var)
    get_trait(tree, "Prior", 1)
end



function Menura_bayesian(reference_data, tree, known_var; trait_drift = trait_drift, trait_diff = trait_diff,
    matrix_drift = covariance_mat_drift, max_iter = convert(Int, 2e6), x0 = nothing, 
    threshold = 0.5, n_particles= 2000)
    #For unknown functions get variables
    priors = JMM_get_priors(tree, known_var)

    # Define function in here - more customisable
    # Defines a version of the simulation function that only takes one input 
    function bayesian_menura_sim(paramaters)
        para = [parameters[1:3], paramaters[4:6], parameters[7:9], paramters[10]]
        res = menura_sim_exper(paramaters..., tree, x0 = x0, trait_drift = trait_drift, trait_diff = trait_diff,
                            matrix_drift = matrix_drift)
        return res[2]
    end
    SimulatedABCRejection(reference_data, bayesian_menura_sim, priors, threshold, n_particles, 
                            max_iter = max_iter, write_progress = true)
end
