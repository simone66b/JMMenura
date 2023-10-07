using GpABC, Distributions, DifferentialEquations, Plots, PosDefManifold, StatsPlots


"""
A function which returns the required priors from a tree

Could make it so you can choose which ones
"""
function JMM_get_priors(tree, known_var)
    get_trait(tree, "Prior", 1)
end



function trait_prior_alpha_all_diff()
    return nothing
end

"""
Function that takes the output of generated priors and returns the interprets it for an alpha value for the traits

"""
function mat_prior_alpha_all_diff(prior_out::Vector{Any})
    n = convert(Int, √length(prior_out))
    return reshape(values, n, n)
end

function prior_inter(;trait_prior_inter = nothing, mat_prior_inter = nothing)
    return nothing
end

# length of priors should be the same. Maybe I need to type this
function extract_priors(trait_priors, mat_priors)
    trait_pri = collect(trait_priors)
    mat_pri = collect(mat_priors)
    collect(Iterators.flatten(append!(trait_pri, mat_pri)))
end


"""
Calculates the distance between two JMMenura simulations. 
"""
function trait_mat_distance(var_num, leaf_num)
    function trait_mat_dist(data1, data2)
        data1 = reshape(data1, var_num, (var_num+1)*leaf_num)
        data2 = reshape(data2, var_num, (var_num+1)*leaf_num)

        traits1 = [data1[:,i] for i in 1:leaf_num]
        traits2 = [data2[:,i] for i in 1:leaf_num]

        mats1 = [data1[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]
        mats2 = [data2[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]

        trait_diff = sum([euclidean(traits1[i], traits2[i]) for i in 1:leaf_num])
        hermi_matrices1, hermi_matrices2 = Hermitian.(mats1), Hermitian.(mats2)
        matrix_diff = sum([PosDefManifold.distance(Fisher, hermi_matrices1[i], hermi_matrices2[i]) for i in 1:leaf_num])
        return trait_diff + matrix_diff
    end
end

"""
Extracts data from a JMMenura simulation into a format usable by GpABC
"""
function get_data(tree)
    traits = reduce(hcat, [tip.data["trait_trace"][end] for tip in getleaves(tree)])
    mats = reduce(hcat, [tip.data["mat_trace"][end] for tip in getleaves(tree)])
    data = [traits..., mats...]
    return reshape(data, length(data), 1)
end 

# Stuff it I'm hard coding an example
function menura_bayesian_trait_alpha(reference_data, tree, trait_known_para, alpha_prior
    ,mat_known_para, trait0, mat0, max_iter, length, error)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)
    
        trait_para = Dict(11 => (trait_known_para..., alpha = parameter.*ones(length))) 

        sim = menura_para_descend!(mat_known_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end

    SimulatedABCRejection(reference_data, bayesian_menura!, [alpha_prior], error, 600, 
    max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(3,6))
end



function menura_bayesian_mat_alpha(reference_data, tree, trait_known_para
    ,mat_known_para, alpha_prior, trait0, mat0, max_iter, length, error)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)
    
        mat_para = Dict(11 => (mat_known_para..., alpha = parameter.* Matrix(1I, length, length))) 

        sim = menura_para_descend!(mat_para, trait_known_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end

    SimulatedABCRejection(reference_data, bayesian_menura!, [alpha_prior], error, 600, 
    max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(3,6))
end

function menura_bayesian_trait_mat_alpha(reference_data, tree, trait_known_para
    ,mat_known_para, alpha_priors, trait0, mat0, max_iter, length, error)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)

        trait_para = Dict(11 => (trait_known_para..., alpha = parameter[1].*ones(length))) 
    
        mat_para = Dict(11 => (mat_known_para..., alpha = parameter[2].* ones(length,length))) 

        sim = menura_para_descend!(mat_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end

    SimulatedABCRejection(reference_data, bayesian_menura!, alpha_priors, error, 600, 
    max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(3,6))
end

"""
Function to perform approximate bayesian computation for using the JMMenura simulation framework. Intakes a series of proposed traits and priors along with reference data.

"""
function menura_bayesian(reference_data, tree, trait_known_para, traits_priors
    ,mat_known_para, mat_priors, trait0, mat0, max_iter, length)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)
    
        trait_para = Dict(11 => (trait_known_para..., alpha = parameter.*ones(length))) 

        sim = menura_para_descend!(mat_known_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end

    SimulatedABCRejection(reference_data, bayesian_menura!, [alpha_prior], 170.0, 600, 
    max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(3,6))
end