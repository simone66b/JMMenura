using DifferentialEquations, Distances, Distributions, JLD2, LinearAlgebra, Phylo, Plots, PyPlot, KissABC

#################################
# Diffusion and Drift functions #
#################################
function diffusion_menura(x0, tspan, p, mat, dt=0.001)
    
    function drift_menura(du, u, p, t)
        alpha = p.alpha;
        mu = p.mu;
        ## du .= alpha .* u would be BM with drift
        du .= alpha .* (mu .- u); ## OU-like
    end; # drift
    
    function diff_menura(du, u, p, t)
        sigma = p.sigma;
        du .= sigma; ## would be OU
        ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
        ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
    end; # diff

    cor1 = cor(mat);
    noise = CorrelatedWienerProcess(cor1, tspan[1],
                                    zeros(dim(cor1)),
                                    zeros(dim(cor1)));
    
    prob = SDEProblem(drift_menura, diff_menura, x0, tspan, p=p, noise=noise);       
    solve(prob, EM(), dt=dt, p=p, adaptive=false);
end; # diffusion

function drift_cov(du, u, p, t) ## drift function for the SDE
    du .= p.a .* t .* p.A;
end # drift

function diffusion_cov(du, u, p, t) ## diffusion function for the SDE
    du .= p.b .* t .* p.B ;
end # diffusion

#######################
# Recursive Functions #
#######################
function recurse_menura!(tree, node, t0 = 0.0)
    if ismissing(node.inbound) ## the root node, to get started
        node.data["trace"]= [x0];
        node.data["timebase"] = [t0];
    else
        ancestor = getancestors(tree, node)[1];
        evol = diffusion_menura(ancestor.data["trace"][end],
                         (getheight(tree, ancestor), getheight(tree, node)),
                         ancestor.data["parameters"], ancestor.data["matrix"]);
       node.data["trace"] = evol.u;
        node.data["timebase"] = evol.t;
        
    end # else
    if !isleaf(tree, node)
        recurse_menura!(tree, node.other[1].inout[2]);
        recurse_menura!(tree, node.other[2].inout[2]);
    end
end # Recurse! 

function recurse_mat!(tree, node)
        
    if ismissing(node.inbound) ## if root
        node.data["matrix"] = node.data["parameters"].mat; ## starting matrix
    else
        ancestor = getancestors(tree, node)[1];
        node.data["matrix"] =
            gen_cov_mat(ancestor.data["matrix"],
                        ancestor.data["parameters"], 
                        (getheight(tree, ancestor),
                         getheight(tree, node)));
    end;
    if !isleaf(tree, node)
        recurse_mat!(tree, node.other[1].inout[2]);
        recurse_mat!(tree, node.other[2].inout[2]);
    end;
end; # Recurse!

############################
# Create Covariance Matrix #
############################

function gen_cov_mat(mat, p, tspan,  u0=zeros(size(mat)), dt = 0.001)
    lowertri = LowerTriangular(mat);
    uppertri = - UpperTriangular(mat);
    skewsymm = lowertri + uppertri;
    W = WienerProcess(0.0, 0.0, 0.0);

    pp = (A=skewsymm, B=skewsymm, a=p.a, b=p.b); ## skew symmetric matrices not necessarily the same.
    prob = SDEProblem(drift_cov, diffusion_cov, u0, tspan, p=pp, noise=W,
                    noise_rate_prototype=zeros(size(mat))); ## setup SDE problem
    sol = solve(prob, EM(), p=pp, dt=dt);
    Omega1 = exp(last(sol.u)); ## get the final matrix

    Omega1 * mat * Omega1'; ## reconstruct P_1
end

##############################
# Trait Prediction Functions #
##############################

function predict_trait_tree(tree)
    #### get the last multivariate trait value in a branch
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

function menura!(tree)
    root = getroot(tree);
    recurse_menura!(tree, root); # do the recursive simulations
    return tree
end # menura!

function menuramat!(tree) ## Only call after putp! 
    root = getroot(tree);
    recurse_mat!(tree, root); # do the recursive simulations
    return tree
end; # menuramat!

function putp!(tree, p1, key)
    for i in eachindex(tree.nodes)
        tree.nodes[i].data[key] = p1;
    end
    tree;
end # putp!

function simulate(parms = parms, mu=mu1, mat=P0, a=a1, b=b1, tree=tree)
    alpha, sigma = parms
    p1 = (alpha=alpha, sigma = sigma, mu=mu, mat = mat, a=a, b=b)
    putp!(tree, p1, "parameters");
    menuramat!(tree);
    menura!(tree);
 (tree, predict_trait_tree(tree))[2];
end