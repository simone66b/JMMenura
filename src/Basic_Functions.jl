using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, PyPlot, KissABC

include("Diffusion_Functions.jl")
include("Evolution_functions.jl")


#######################
# Recursive Functions #
#######################
function recurse_menura!(tree, node, t0 , x0)
    if ismissing(node.inbound) ## the root node, to get started
        node.data["trace"]= [x0]
        node.data["timebase"] = [t0]
    else
        ancestor = getancestors(tree, node)[1]
        evol = trait_diffusion(ancestor.data["trace"][end],
                         (getheight(tree, ancestor), getheight(tree, node)),
                         ancestor.data["parameters"], ancestor.data["matrix"])
        node.data["trace"] = evol.u
        node.data["timebase"] = evol.t
        
    end # else
    if !isleaf(tree, node)
        recurse_menura!(tree, node.other[1].inout[2], t0 , x0)
        recurse_menura!(tree, node.other[2].inout[2], t0 , x0)
    end
end # Recurse! 

function recurse_mat!(tree, node)
        
    if ismissing(node.inbound) ## if root
        node.data["matrix"] = node.data["parameters"].mat ## starting matrix
    else
        ancestor = getancestors(tree, node)[1]
        node.data["matrix"] =
            gen_cov_mat(ancestor.data["matrix"],
                        ancestor.data["parameters"], 
                        (getheight(tree, ancestor),
                         getheight(tree, node)))
    end;
    if !isleaf(tree, node)
        recurse_mat!(tree, node.other[1].inout[2])
        recurse_mat!(tree, node.other[2].inout[2])
    end;
end; # Recurse!

##############################
# Trait Prediction Functions #
##############################

function predict_trait_tree(tree)
    """
    get the last multivariate trait value in a branch
    """
    testtips = getleaves(tree)
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

function menura!(tree, para_len)
    # println("length = $para_len")
    x0 = repeat([0.0], para_len)
    root = getroot(tree);
    recurse_menura!(tree, root, 0.0, x0) # do the recursive simulations
    return tree
end # menura!

function menuramat!(tree) ## Only call after putp! 
    root = getroot(tree)
    recurse_mat!(tree, root) # do the recursive simulations
    return tree
end; # menuramat!

function putp!(tree, p1, key)
    for i in eachindex(tree.nodes)
        tree.nodes[i].data[key] = p1
    end
    tree;
end # putp!

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

function menura_sim(alpha, sigma, mu, mat, a, b, tree)
    """
    Hello
    """
    #menura_errors(alpha, sigma, mu, mat)

    para_len = length(alpha)

    p1 = (alpha=alpha, sigma=sigma, mu=mu, mat=mat, a=a, b=b)
    putp!(tree, p1, "parameters")
    menuramat!(tree)
    menura!(tree, para_len)
 (tree, predict_trait_tree(tree))[2]
end