

# Simulation functtio
function JMMABC(reference_data, tree, JMMpara::JMMABCparameters, trait0, mat0, threshold, n_particles; max_iter = 50*n_particles, t0 = 0.0, each = false, 
    dt = 0.001, distance_function = trait_mat_distance(JMMpara.size,nleaves(tree)), summary_function = get_data, verbose = true) # OTHER things needed for simulation
    # get priors

    bayesian_menura! = create_bayesian_sim(tree, JMMpara, trait0, mat0, t0 = t0, each = each, dt = dt, summary_function = summary_function, verbose = verbose)

    while 
end
