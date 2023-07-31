using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, PyPlot, KissABC
#######################
# Recursive Functions #
#######################

"""
Recursively iterates over the tree and evolves each node.
"""
function recurse_menura_each!(tree, node, t0 , x0, trait_drift, trait_diff, matrix_drift, matrix_func; 
                                small_dt_scale = 10)
    if ismissing(node.inbound) ## the root node, to get started
        node.data["matrixes"] = [node.data["parameters"].mat] ## starting matrix

        node.data["trace"]= [x0]
        node.data["timebase"] = [t0]
    else   
        ancestor = getancestors(tree, node)[1]

        if haskey(node.data, "known_G_mat")
            evol_matrix = node.data["known_G_mat"]
         else 
             evol_matrix = ancestor.data["matrixes"][end]
         end             
        node.data["matrixes"] =
            matrix_func(evol_matrix,
                        ancestor.data["parameters"], 
                        (getheight(tree, ancestor),
                         getheight(tree, node)), matrix_drift)

        evol = trait_diffusion_each(ancestor.data["trace"][end],
                         (getheight(tree, ancestor), getheight(tree, node)),
                         ancestor.data["parameters"], 
                         node.data["matrixes"], trait_drift, trait_diff, small_dt_scale = small_dt_scale)
        node.data["trace"] = evol.u
        node.data["timebase"] = evol.t
        
    end
    if !isleaf(tree, node)
        recurse_menura_each!(tree, node.other[1].inout[2], t0 , x0, trait_drift, trait_diff, matrix_drift
                                , matrix_func, small_dt_scale = small_dt_scale)
        recurse_menura_each!(tree, node.other[2].inout[2], t0 , x0, trait_drift, trait_diff, matrix_drift
                                , matrix_func, small_dt_scale = small_dt_scale)
    end
end


##############################
# Trait Prediction Functions #
##############################

# """
#     predict_trait_tree(tree)

# Returns the last multivariate trait value from all branches
# """
# function predict_trait_tree(tree)
#     testtips = getleaves(tree);
#     res = Array{Vector{Float64}}(undef, length(testtips))
#     tipnames = Array{String}(undef, length(testtips))
#     tiptimes = Vector{Float64}()
    
#     for i in eachindex(testtips) ## could maybe use heightstoroot() for this computation
#         res[i] = testtips[i].data["trace"][end]
#         tipnames[i] = testtips[i].name;
#         push!(tiptimes, getheight(tree, testtips[i]))
#     end;
#     collect(Iterators.flatten(res))
# end # predictTraitTree

#################################
# Simulation Handling Functions #
#################################

"""
Now that i think about it this function is possibly redundant.
"""
function menura_each!(tree, x0, trait_drift, trait_diff, matrix_drift, matrix_func; small_dt_scale = 10)
    root = getroot(tree)
    recurse_menura_each!(tree, root, 0.0, x0, trait_drift, trait_diff, matrix_drift, matrix_func,
                            small_dt_scale = small_dt_scale )
end # menura_each!

# """
# Makes sure input parameters are acceptable (Might need to refine more).
# """
# function menura_errors(alpha, sigma, mu, mat)
#     if size(mat)[1] != size(mat)[2]
#         error("Covariance matrix must be square matrix")
#     end

#     if length(alpha) != length(sigma)
#         error("Alpha and Sigma must be the same length")
#     elseif length(mu) != length(sigma)
#         error("Mu and Sigma must be the same length")
#     elseif length(mu) != length(alpha)
#         error("Alpha and Sigma must be the same length")
#     elseif length(mu) != size(mat)[1]
#         error("Mu and Covariance matrix must be the same size")
#     elseif length(alpha) != size(mat)[1]
#         error("Mu and Covariance matrix must be the same size")
#     elseif length(sigma) != size(mat)[1]
#         error("Mu and Covariance matrix must be the same size")
#     end
# end

"""
    menura_sim(alpha, sigma, mu, mat, a, b, tree; x0 = nothing, trait_drift = trait_drift, trait_diff = trait_diff,
                        matrix_drift = covariance_mat_drift)

Handles simulating evolution over a phylogenetic tree.

Returns the tree, along with the trait values on each of the leaves of the tree.

# Arguments 
- 
"""
function menura_sim_each(alpha, sigma, mu, cov_mat, tree, matrix_func; a = nothing, b = nothing, x0 = nothing, mat_alpha = nothing
                            , mat_sigma = nothing, mat_mu = nothing, trait_drift = trait_drift, trait_diff = trait_diff, 
                            matrix_drift = covariance_mat_drift, small_dt_scale = 10)
    # menura_errors(alpha, sigma, mu, cov_mat)

    para_len = length(alpha)
    if isnothing(x0)
        x0 = repeat([0.0], para_len)
    end

    p1 = (alpha=alpha, sigma=sigma, mu=mu, mat=cov_mat, a=a, b=b, mat_alpha = mat_alpha, mat_sigma = mat_sigma, mat_mu = mat_mu)
    putp!(tree, p1, "parameters")
    menura_each!(tree, x0, trait_drift, trait_diff, matrix_drift, matrix_func, small_dt_scale = small_dt_scale )
    # (tree, predict_trait_tree(tree))
    (tree, reduce(hcat, [tip.data["trace"][end] for tip in getleaves(tree)]))
end

# function menura_sim_mat_isospectral_each(alpha, sigma, mu, cov_mat, a, b, tree; x0 = nothing, trait_drift = trait_drift, trait_diff = trait_diff,
#     matrix_drift = covariance_mat_drift)
#     menura_sim_exper(alpha, sigma, mu, cov_mat, tree, gen_cov_mat, a = a, b = b, x0 = x0, trait_drift = trait_drift, trait_diff = trait_diff, 
#                             matrix_drift = covariance_mat_drift)
# end

function menura_sim_mat_OU_each(alpha, sigma, mu, cov_mat, mat_alpha, mat_sigma, mat_mu, tree; 
    x0 = nothing, trait_drift = trait_drift, trait_diff = trait_diff,
    matrix_drift = covariance_mat_drift, small_dt_scale = 10)
    menura_sim_each(alpha, sigma, mu, cov_mat, tree, OUmatrix_each, mat_alpha = mat_alpha, mat_sigma = mat_sigma, mat_mu = mat_mu, x0 = x0, trait_drift = trait_drift, trait_diff = trait_diff, 
                            matrix_drift = matrix_OU_drift, small_dt_scale = small_dt_scale)
end