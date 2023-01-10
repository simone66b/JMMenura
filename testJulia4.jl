using DifferentialEquations, Phylo, Plots, Distributions, Distances, JLD2, LinearAlgebra, GpABC;
pyplot();

    #####################################################################################
    #####################################################################################
    
    function menura!(tree)

        function diffusion(x0, tspan, p, mat, dt=0.001)
            function drift(du, u, p, t)
                alpha = p.alpha;
                mu = p.mu;
                ## du .= alpha .* u would be BM with drift
                du .= alpha .* (mu .- u); ## OU-like
            end; # drift
            
            function diff(du, u, p, t)
                sigma = p.sigma;
                du .= sigma; ## would be OU
                ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
                ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
            end; # diff
            
            cor1 = cor(mat);
            noise = CorrelatedWienerProcess(cor1, tspan[1],
                                            zeros(dim(cor1)),
                                            zeros(dim(cor1)));
            
            prob = SDEProblem(drift, diff, x0, tspan, p=p, noise=noise);       
            solve(prob, EM(), dt=dt, p=p, adaptive=false);
        end; # diffusion

        ################################################################################
        
        function Recurse!(tree, node, t0 = 0.0)
            if ismissing(node.inbound) ## the root node, to get started
                node.data["trace"]= [x0];
                node.data["timebase"] = [t0];
            else
                ancestor = getancestors(tree, node)[1];
                evol = diffusion(ancestor.data["trace"][end],
                                 (getheight(tree, ancestor), getheight(tree, node)),
                                 ancestor.data["parameters"], ancestor.data["matrix"]);
               node.data["trace"] = evol.u;
                node.data["timebase"] = evol.t;
                
            end # else
            if !isleaf(tree, node)
                Recurse!(tree, node.other[1].inout[2]);
                Recurse!(tree, node.other[2].inout[2]);
            end
        end # Recurse!  
        root = getroot(tree);
        Recurse!(tree, root); # do the recursive simulations
        tree;
    end # menura!
    
    ###############################################################################
    ################################################################################
    
    function menuramat!(tree) ## Only call after putp! 
        function Recurse!(tree, node)
            
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
                Recurse!(tree, node.other[1].inout[2]);
                Recurse!(tree, node.other[2].inout[2]);
            end;
        end; # Recurse!
        
        root = getroot(tree);
        Recurse!(tree, root); # do the recursive simulations
        tree;
    end; # menuramat!
    
    ############################################################################33
    #############################################################################
    
    function predictTraitTree(tree)
        #### get the last multivariate trait value in a branch
        testtips = getleaves(tree);
        res = Array{Vector{Float64}}(undef, length(testtips)); 
        tipnames = Array{String}(undef, length(testtips));
        tiptimes = Vector{Float64}();
        ##  finaltraitvals = Dict();
        ## tiptimesdict = Dict();
        
        for i in 1:length(testtips) ## could maybe use heightstoroot() for this computation
            res[i] = testtips[i].data["trace"][end];
            tipnames[i] = testtips[i].name;
            push!(tiptimes, getheight(tree, testtips[i]));
        end;
        collect(Iterators.flatten(res))
    end # predictTraitTree
    
    ##########################################################################33
    ###########################################################################3
    
    function putp!(tree, p1, key)
        for i in 1:length(tree.nodes)
            tree.nodes[i].data[key] = p1;
        end
        tree;
    end # putp!
    
    ####################################################################33
    ######################################################################3
function gen_cov_mat(mat, p, tspan,  u0=zeros(size(mat)), dt = 0.001)
    function drift(du, u, p, t) ## drift function for the SDE
        du .= p.a .* t .* p.A;
    end # drift
    
    function diffusion(du, u, p, t) ## diffusion function for the SDE
        du .= p.b .* t .* p.B ;
    end # diffusion
    
    lowertri = LowerTriangular(mat);
    uppertri = - UpperTriangular(mat);
    skewsymm = lowertri + uppertri;
    W = WienerProcess(0.0, 0.0, 0.0);
    
    pp = (A=skewsymm, B=skewsymm, a=p.a, b=p.b); ## skew symmetric matrices not necessarily the same.
    prob = SDEProblem(drift, diffusion, u0, tspan, p=pp, noise=W,
                      noise_rate_prototype=zeros(size(mat))); ## setup SDE problem
    sol = solve(prob, EM(), p=pp, dt=dt);
    Omega1 = exp(last(sol.u)); ## get the final matrix

    Omega1 * mat * Omega1'; ## reconstruct P_1
end # gen_cov_mat
    ######################################################################3
    ##################################################################3n

    a1=1.0;
    b1= 1.0;
x0 =  repeat([0.0], 8);
tree1 = Ultrametric(20);
    tree = rand(tree1); 
       time_tot = 1.0;
    tspan = (0.0, time_tot);

P0 = [0.329 0.094 -0.083 -0.089 0.293 0.079 0.208 0.268;
0.094 0.449 0.349 0.24 0.071 0.075 0.03 0.009;
-0.083 0.349 1.426 0.487 -0.371 -0.098 -0.053 -0.172;
-0.089 0.24 0.487 0.546 -0.168 0.017 -0.051 -0.081;
0.293 0.071 -0.371 -0.168 1.441 1.008 0.904 0.945;
0.079 0.075 -0.098 0.017 1.008 1.087 0.731 0.78;
0.208 0.03 -0.053 -0.051 0.904 0.731 0.809 0.783;
0.268 0.009 -0.172 -0.081 0.945 0.78 0.783 0.949];

alpha1 = repeat([1.0], 8);
mu1 = repeat([0.0], 8); ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], 8);

function simulate(tree, alpha=alpha1, sigma=sigma1, mu=mu1, mat=P0, a=a1, b=b1)
    p1 = (alpha=alpha, mu=mu, sigma=sigma, mat=P0, a=a, b=b);
    putp!(tree, p1, "parameters");
    menuramat!(tree);
    menura!(tree);
(tree, predictTraitTree(tree));
end; # simulate

exampledat = simulate(tree, alpha1, sigma1, mu1);

# priordists = [Truncated(Normal(0,3), 0, Inf)]
# priors = Factored(repeat(priordists, 12)...,)

# function cost((alpha, mu, sigma))
#     x=simulate((alpha, mu, sigma))[2]
#     y=exampledat[2]
#     euclidean(x, y)
# end #cost

# approx_density = ApproxKernelizedPosterior(priors,cost,0.005)
# ## res = sample(approx_density, AIS(25), 1000, ntransitions=100, discard_initial = 10)   

## save_object("ABCResults.jld2", res)


########################################################################
#########################################################################

 plot(tree,
    size = (400, 800),
    ## markersize = 20, 
      series_annotations = text.(1:nnodes(tree), 15, :center, :center,
                                 :white),  linewidth=5, showtips=false)
current()

testnodes = getnodes(exampledat[1]);

plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0),
     legend=nothing,
     reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 8, length(u1)));
    myt = i.data["timebase"];
    plot!(myt, uu1[:, 1], uu1[:,3]);
end; # for
current()

plot(xlim=(0.0,1.0), ylim= (-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 8, length(u1)));
    myt = i.data["timebase"];
    plot!(myt, uu1[:, 3]);
end; # for
current()

plot(ylim=(-2.0, 2.0), xlim= (-2.0, 2.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 8, length(u1)));
    plot!(uu1[:,1], uu1[:, 3]);
end;
current()

plot(xlim=(-2.0,2.0), ylim= (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing,
     reuse=false)
for i in testnodes
    u1= i.data["trace"];
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 8, length(u1)));
    myt = i.data["timebase"];
    plot!(uu1[:,3], uu1[:, 1], uu1[:,2]);
end;
current()



######################################################################333
#########################################################################
##########################################################################
