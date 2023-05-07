using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, PyPlot, KissABC

include("Diffusion_Functions.jl")



#######################
# Recursive Functions #
#######################
function recurse_menura!(tree, node, t0 , x0)
    if ismissing(node.inbound) ## the root node, to get started
        node.data["matrix"] = node.data["parameters"].mat; ## starting matrix

        node.data["trace"]= [x0];
        node.data["timebase"] = [t0];
    else
        ancestor = getancestors(tree, node)[1];
        node.data["matrix"] =
            gen_cov_mat(ancestor.data["matrix"],
                        ancestor.data["parameters"], 
                        (getheight(tree, ancestor),
                         getheight(tree, node)));

        evol = diffusion_menura(ancestor.data["trace"][end],
                         (getheight(tree, ancestor), getheight(tree, node)),
                         ancestor.data["parameters"], ancestor.data["matrix"]);
        node.data["trace"] = evol.u
        node.data["timebase"] = evol.t
        
    end # else
    if !isleaf(tree, node)
        recurse_menura!(tree, node.other[1].inout[2], t0 , x0);
        recurse_menura!(tree, node.other[2].inout[2], t0 , x0);
    end
end # Recurse! 


##############################
# Matrix and trait evolution #
##############################

function diffusion_menura(x0, tspan, p, mat, dt=0.001)
    cor1 = cor(mat);
    noise = CorrelatedWienerProcess(cor1, tspan[1],
                                    zeros(dim(cor1)),
                                    zeros(dim(cor1))); #Make a length function
    
    prob = SDEProblem(trait_drift, trait_diff, x0, tspan, p=p, noise=noise);       
    solve(prob, EM(), dt=dt, p=p, adaptive=false);
end; # diffusion

function gen_cov_mat(mat, p, tspan,  u0=zeros(size(mat)), dt = 0.001)
    lowertri = LowerTriangular(mat);
    uppertri = - UpperTriangular(mat);
    skewsymm = lowertri + uppertri;
    W = WienerProcess(0.0, 0.0, 0.0);

    pp = (A=skewsymm, B=skewsymm, a=p.a, b=p.b); ## skew symmetric matrices not necessarily the same.
    prob = SDEProblem(covariance_mat_drift, covariance_mat_diffusion, u0, tspan, p=pp, noise=W,
                    noise_rate_prototype=zeros(size(mat))); ## setup SDE problem
    sol = solve(prob, EM(), p=pp, dt=dt);
    Omega1 = exp(last(sol.u)); ## get the final matrix

    Omega1 * mat * Omega1'; ## reconstruct P_1
end

##############################
# Trait Prediction Functions #
##############################

function predict_trait_tree(tree)
    """
    Get the last multivariate trait value from all branches
    """
    testtips = getleaves(tree);
    res = Array{Vector{Float64}}(undef, length(testtips)); 
    tipnames = Array{String}(undef, length(testtips));
    tiptimes = Vector{Float64}();
    
    for i in eachindex(testtips) ## could maybe use heightstoroot() for this computation
        res[i] = testtips[i].data["trace"][end];
        tipnames[i] = testtips[i].name;
        push!(tiptimes, getheight(tree, testtips[i]));
    end;
    collect(Iterators.flatten(res))
end # predictTraitTree

#################################
# Simulation Handling Functions #
#################################

function menura!(tree, para_len, x0)
    root = getroot(tree);
    recurse_menura!(tree, root, 0.0, x0); # do the recursive simulations
    return tree
end # menura!

function putp!(tree, p1, key)
    for i in eachindex(tree.nodes)
        tree.nodes[i].data[key] = p1;
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

function menura_sim_exper(alpha, sigma, mu, mat, a, b, tree; x0 = nothing)
    """
    Hello
    """
    #menura_errors(alpha, sigma, mu, mat)

    para_len = length(alpha)
    if isnothing(x0)
        x0 = repeat([0.0], para_len)
    end

    p1 = (alpha=alpha, sigma=sigma, mu=mu, mat=mat, a=a, b=b)
    putp!(tree, p1, "parameters");
    menura!(tree, para_len, x0);
    (tree, predict_trait_tree(tree))[2];
end