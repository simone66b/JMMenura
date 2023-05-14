using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, PyPlot, KissABC

include("Diffusion_Functions.jl")
include("Evolution_functions.jl")
include("Tree_Modifying_Functions.jl")


#######################
# Recursive Functions #
#######################
function recurse_menura!(tree, node, t0 , x0, trait_drift, trait_diff, matrix_drift)
    if ismissing(node.inbound) ## the root node, to get started
        node.data["matrix"] = node.data["parameters"].mat ## starting matrix

        node.data["trace"]= [x0]
        node.data["timebase"] = [t0]
    else
        ancestor = getancestors(tree, node)[1]
        node.data["matrix"] =
            gen_cov_mat(ancestor.data["matrix"],
                        ancestor.data["parameters"], 
                        (getheight(tree, ancestor),
                         getheight(tree, node)), matrix_drift)

        evol = trait_diffusion(ancestor.data["trace"][end],
                         (getheight(tree, ancestor), getheight(tree, node)),
                         ancestor.data["parameters"], ancestor.data["matrix"], trait_drift, trait_diff)
        node.data["trace"] = evol.u
        node.data["timebase"] = evol.t
        
    end # else
    if !isleaf(tree, node)
        recurse_menura!(tree, node.other[1].inout[2], t0 , x0, trait_drift, trait_diff, matrix_drift)
        recurse_menura!(tree, node.other[2].inout[2], t0 , x0, trait_drift, trait_diff, matrix_drift)
    end
end # Recurse! 


##############################
# Trait Prediction Functions #
##############################

function predict_trait_tree(tree)
    """
    Get the last multivariate trait value from all branches
    """
    testtips = getleaves(tree);
    res = Array{Vector{Float64}}(undef, length(testtips))
    tipnames = Array{String}(undef, length(testtips))
    tiptimes = Vector{Float64}()
    
    for i in eachindex(testtips) ## could maybe use heightstoroot() for this computation
        res[i] = testtips[i].data["trace"][end]
        tipnames[i] = testtips[i].name;
        push!(tiptimes, getheight(tree, testtips[i]))
    end;
    collect(Iterators.flatten(res))
end # predictTraitTree

#################################
# Simulation Handling Functions #
#################################

function menura!(tree, x0, trait_drift, trait_diff, matrix_drift)
    root = getroot(tree)
    recurse_menura!(tree, root, 0.0, x0, trait_drift, trait_diff, matrix_drift)
end # menura!

function menura_errors(alpha, sigma, mu, mat)
    if size(mat)[1] != size(mat)[2]
        error("Covariance matrix must be square matrix")
    end

    if length(alpha) != length(sigma)
        error("Alpha and Sigma must be the same length")
    elseif length(mu) != length(sigma)
        error("Mu and Sigma must be the same length")
    elseif length(mu) != length(alpha)
        error("Alpha and Sigma must be the same length")
    elseif length(mu) != size(mat)[1]
        error("Mu and Covariance matrix must be the same size")
    elseif length(alpha) != size(mat)[1]
        error("Mu and Covariance matrix must be the same size")
    elseif length(sigma) != size(mat)[1]
        error("Mu and Covariance matrix must be the same size")
    end
end

function menura_sim_exper(alpha, sigma, mu, mat, a, b, tree; x0 = nothing, trait_drift = trait_drift, trait_diff = trait_diff,
                            matrix_drift = covariance_mat_drift)
    """
    Hello
    """
    menura_errors(alpha, sigma, mu, mat)

    para_len = length(alpha)
    if isnothing(x0)
        x0 = repeat([0.0], para_len)
    end

    p1 = (alpha=alpha, sigma=sigma, mu=mu, mat=mat, a=a, b=b)
    putp!(tree, p1, "parameters")
    menura!(tree, x0, trait_drift, trait_diff, matrix_drift)
    (tree, predict_trait_tree(tree))
end