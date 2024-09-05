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
function trait_mat_distance(var_num, leaf_num; err_thres = 10^-14)
    function trait_mat_dist(data1, data2)
        data1 = reshape(data1, var_num, (var_num+1)*leaf_num)
        data2 = reshape(data2, var_num, (var_num+1)*leaf_num)

        traits1 = [data1[:,i] for i in 1:leaf_num]
        traits2 = [data2[:,i] for i in 1:leaf_num]

        mats1 = [data1[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]
        mats2 = [data2[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]

        trait_diff =  mean((traits1[i][j] - traits2[i][j])^2 for i in 1:leaf_num for j in 1:var_num)^0.5
        hermi_matrices1, hermi_matrices2 = Hermitian.(mats1), Hermitian.(mats2)

        # Stupid floating points
        hermi_matrices1 = correct_mat.(hermi_matrices1, err_thres)
        hermi_matrices2 = correct_mat.(hermi_matrices2, err_thres)

        # matrix_diff = mean([PosDefManifold.distance(Fisher, hermi_matrices1[i], hermi_matrices2[i]) for i in 1:leaf_num])
        matrix_diff = mean([sqrt(sum(log.(max.(eigvals(hermi_matrices1[i], hermi_matrices2[i]), 0)).^2)) for i in 1:leaf_num]) # Fisher Rao metric
        return trait_diff + matrix_diff
    end
end

function trait_mat_distance_scaled(var_num, leaf_num, mat_scale = 1, trait_scale = 1; err_thres = 10^-14)
    function trait_mat_dist(data1, data2)
        data1 = reshape(data1, var_num, (var_num+1)*leaf_num)
        data2 = reshape(data2, var_num, (var_num+1)*leaf_num)

        traits1 = [data1[:,i] for i in 1:leaf_num]
        traits2 = [data2[:,i] for i in 1:leaf_num]

        mats1 = [data1[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]
        mats2 = [data2[:,(leaf_num + var_num*(i-1)+1):(leaf_num + var_num*(i))] for i in 1:leaf_num]

        trait_diff =  mean((traits1[i][j] - traits2[i][j])^2 for i in 1:leaf_num for j in 1:var_num)^0.5
        hermi_matrices1, hermi_matrices2 = Hermitian.(mats1), Hermitian.(mats2)

        # Stupid floating points
        hermi_matrices1 = correct_mat.(hermi_matrices1, err_thres)
        hermi_matrices2 = correct_mat.(hermi_matrices2, err_thres)

        # matrix_diff = mean([PosDefManifold.distance(Fisher, hermi_matrices1[i], hermi_matrices2[i]) for i in 1:leaf_num])
        matrix_diff = mean([sqrt(sum(log.(max.(eigvals(hermi_matrices1[i], hermi_matrices2[i]), 0)).^2)) for i in 1:leaf_num]) # Fisher Rao metric
        return trait_scale*trait_diff + mat_scale*matrix_diff
    end
end

function mat_distance(var_num, leaf_num; err_thres = 10^-14)
    function mat_dist(data1, data2)
        data1 = reshape(data1, var_num, (var_num)*leaf_num)
        data2 = reshape(data2, var_num, (var_num)*leaf_num)

        traits1 = [data1[:,i] for i in 1:leaf_num]
        traits2 = [data2[:,i] for i in 1:leaf_num]

        mats1 = [data1[:,(var_num*(i-1)+1):(var_num*(i))] for i in 1:leaf_num]
        mats2 = [data2[:,(var_num*(i-1)+1):(var_num*(i))] for i in 1:leaf_num]

        hermi_matrices1, hermi_matrices2 = Hermitian.(mats1), Hermitian.(mats2)

        # Stupid floating points
        hermi_matrices1 = correct_mat.(hermi_matrices1, err_thres)
        hermi_matrices2 = correct_mat.(hermi_matrices2, err_thres)

        # matrix_diff = mean([PosDefManifold.distance(Fisher, hermi_matrices1[i], hermi_matrices2[i]) for i in 1:leaf_num])
        matrix_diff = mean([sqrt(sum(log.(max.(eigvals(hermi_matrices1[i], hermi_matrices2[i]), 0)).^2)) for i in 1:leaf_num]) # Fisher Rao metric
        return matrix_diff
    end
end

function trait_distance(var_num, leaf_num; err_thres = 10^-14)
    function trait_dist(data1, data2)
        data1 = reshape(data1, var_num, leaf_num)
        data2 = reshape(data2, var_num, leaf_num)

        traits1 = [data1[:,i] for i in 1:leaf_num]
        traits2 = [data2[:,i] for i in 1:leaf_num]


        trait_diff =  mean((traits1[i][j] - traits2[i][j])^2 for i in 1:leaf_num for j in 1:var_num)^0.5

        return trait_diff
    end
end

"""
Attempt to fix floating point errors
"""
# function correct_mat(mat)
#     println("Eigen: ", eigen(mat).values, "\n")
#     err = eigen(mat).values[1]
#     println("Error: ", err,"\n", isposdef(mat), "\n")
#     if err < 0
#         mat_err_mat = Matrix((err)I, size(mat)...)
#         println("Corrected eigen: ", eigen(mat - mat_err_mat).values, "\n")
#         mat = Hermitian(mat - mat_err_mat)
#         println("Hermitian corrected eigen: ", eigen(mat).values, "\n", isposdef(mat), "\n")
#     end
#     return mat
# end

function correct_mat(mat, err_thres)
    # println("Eigen: ", eigen(mat).values, "\n")
    err = eigen(mat).values[1]
    # println("Error: ", err,"\n", isposdef(mat), "\n")
    if err < err_thres
        mat_err_mat = Matrix((min(-10^-14, err))I, size(mat)...)
        # println("Corrected eigen: ", eigen(mat - mat_err_mat).values, "\n")
        mat = Hermitian(mat - mat_err_mat)
        # println("Hermitian corrected eigen: ", eigen(mat).values, "\n", isposdef(mat), "\n")
    end
    return mat
end

"""
Extracts data from a JMMenura simulation into a format usable by GpABC
"""
function get_data(sim_data)
    tree = sim_data[1]
    traits = reduce(hcat, [tip.data["trait_trace"][end] for tip in getleaves(tree)])
    mats = reduce(hcat, [tip.data["mat_trace"][end] for tip in getleaves(tree)])
    data = [traits..., mats...]
    return sim_data[2] ? reshape(data, length(data), 1) : Array{Float64}(undef, 0, 0)
end 

function get_data_no_trait(sim_data)
    tree = sim_data[1]
    mats = reduce(hcat, [tip.data["mat_trace"][end] for tip in getleaves(tree)])
    return sim_data[2] ? reshape(mats, length(mats), 1) : Array{Float64}(undef, 0, 0)
end 

function get_data_no_mat(sim_data)
    tree = sim_data[1]
    traits = reduce(hcat, [tip.data["trait_trace"][end] for tip in getleaves(tree)])
    return sim_data[2] ? reshape(traits, length(traits), 1) : Array{Float64}(undef, 0, 0)
end 


# Stuff it I'm hard coding an example
function menura_bayesian_trait_alpha(reference_data, tree, trait_known_para, alpha_prior
    ,mat_known_para, trait0, mat0, max_iter, length, error)

    # Pull out priors and place into ordered vector 
    function bayesian_menura!(parameter)
    
        trait_para = Dict(11 => (trait_known_para..., alpha = parameter.*ones(length))) 

        sim = menura_parameter_descend!(mat_known_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

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

        sim = menura_parameter_descend!(mat_para, trait_known_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

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

        sim = menura_parameter_descend!(mat_para, trait_para, tree, trait_evol(), mat_evol(), 0.0, trait0, mat0, false)

        return get_data(sim[1])
    end

    SimulatedABCRejection(reference_data, bayesian_menura!, alpha_priors, error, 600, 
    max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(3,6))
end


"""
Function to perform approximate bayesian computation for using the JMMenura simulation framework. Intakes a series of proposed traits and priors along with reference data.

"""
function menura_bayesian(reference_data, tree, JMMpara::JMMABCparameters, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
    dt = 0.001, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
    
    preallocate_tree!(tree, dt, JMMpara.size)

    bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

    SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
    max_iter = max_iter, write_progress = true, distance_function = distance_function)
end

function test_threshold(reference_data, tree, JMMpara::JMMABCparameters, trait0, mat0, n_particles; 
    t0 = 0.0, each = false, dt = 0.001, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))

    preallocate_tree!(tree, dt, JMMpara.size)

    bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

    thresholds = zeros(n_particles)

    for i in ProgressBar(1:n_particles)
        para = rand.(get_priors(JMMpara))
        sim = bayesian_menura!(para)
        if length(sim) != 0
            dist = distance_function(sim, reference_data)
            thresholds[i] = dist
        end
    end
    return thresholds
end


# function menura_bayesian(reference_data, tree, JMMpara::JMMABCAlphaEqualConstant, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end



# function test_threshold(reference_data, tree, JMMpara::JMMABCAlphaEqualConstant, trait0, mat0, n_particles
#                         ; t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end


# function menura_bayesian(reference_data, tree, JMMpara::JMMABCAlphaDifferentConstant, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     # Pull out priors and place into ordered vector 
#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end



# function test_threshold(reference_data, tree, JMMpara::JMMABCAlphaDifferentConstant, trait0, mat0, n_particles; 
#                         t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end


# function menura_bayesian(reference_data, tree, JMMpara::JMMABCAlphaConstantEqual, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     # Pull out priors and place into ordered vector 
#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end

# function test_threshold(reference_data, tree, JMMpara::JMMABCAlphaConstantEqual, trait0, mat0, n_particles;
#                          t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end

# function menura_bayesian(reference_data, tree, JMMpara::JMMABCIsospectralAlpha, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     # Pull out priors and place into ordered vector 
#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end



# function test_threshold(reference_data, tree, JMMpara::JMMABCIsospectralAlpha, trait0, mat0, n_particles; 
#                         t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end


# function menura_bayesian(reference_data, tree, JMMpara::JMMABCIsospectralAlphaAB, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     # Pull out priors and place into ordered vector 
#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end



# function test_threshold(reference_data, tree, JMMpara::JMMABCIsospectralAlphaAB, trait0, mat0, n_particles; 
#                         t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end

# function menura_bayesian(reference_data, tree, JMMpara::JMMABCBrownian, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
#     dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     # Pull out priors and place into ordered vector
#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     SimulatedABCSMC(reference_data, bayesian_menura!, get_priors(JMMpara), [threshold], n_particles, 
#     max_iter = max_iter, write_progress = true, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)))
# end



# function test_threshold(reference_data, tree, JMMpara::JMMABCBrownian, trait0, mat0, n_particles; 
#                         t0 = 0.0, each = false, dt = 0.001)

#     preallocate_tree!(tree, dt, JMMpara.size)

#     bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt)

#     thresholds = zeros(n_particles)

#     for i in ProgressBar(1:n_particles)
#         para = rand.(get_priors(JMMpara))
#         sim = bayesian_menura!(para)
#         dist = trait_mat_distance(JMMpara.size,nleaves(tree))(sim, reference_data)
#         thresholds[i] = dist
#     end
#     return thresholds
# end