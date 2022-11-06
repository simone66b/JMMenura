using DifferentialEquations, Phylo, Plots, Distributions, Distances, KissABC, JLD2, LinearAlgebra; pyplot();

    #####################################################################################
    #####################################################################################
    
    function menura!(tree)

        function diffusion(x0, tspan, p, mat, dt=0.001)
            function drift(du, u, p, t)
                alpha = p.alpha;
                mu = p.mu;
                ## du .= alpha .* u would be BM with drift
                du .= alpha .* (mu .- u); ## OU-like
            end # drift
            
            function diff(du, u, p, t)
                sigma = p.sigma;
                du .= sigma; ## would be OU
                ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
                ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
            end # diff
            
            cor1 = cor(mat);
            noise = CorrelatedWienerProcess(cor1, tspan[1],
                                            zeros(dim(cor1)),
                                            zeros(dim(cor1)));
            
            prob = SDEProblem(drift, diff, x0, tspan, p=p, noise=noise);       
            solve(prob, EM(), dt=dt, p=p, adaptive=false);
        end # diffusion

        ################################################################################
        
        function Recurse!(tree, node, t0 = 0.0)
            if ismissing(node.inbound) ## the root node, to get started
                node.data["trace"]= [x0];
                node.data["timebase"] = [t0];
            else
                ancestor = getancestors(tree, node)[1]
                evol = diffusion(ancestor.data["trace"][end],
                                 (getheight(tree, ancestor), getheight(tree, node)),
                                 ancestor.data["parameters"], ancestor.data["matrix"]);
               node.data["trace"] = evol.u
                node.data["timebase"] = evol.t
                
            end # else
            if !isleaf(tree, node)
                Recurse!(tree, node.other[1].inout[2]);
                Recurse!(tree, node.other[2].inout[2]);
            end
        end # Recurse!  
        root = getroot(tree)
        Recurse!(tree, root) # do the recursive simulations
        tree
    end # menura!
    
    ###############################################################################
    ################################################################################
    
    function menuramat!(tree) ## Only call after putp! 
        function Recurse!(tree, node)
            
            if ismissing(node.inbound) ## if root
                node.data["matrix"] = node.data["parameters"].mat ## starting matrix
            else
                ancestor = getancestors(tree, node)[1]
                node.data["matrix"] =
                    gen_cov_mat(ancestor.data["matrix"],
                                ancestor.data["parameters"], 
                                (getheight(tree, ancestor),
                                 getheight(tree, node)));
            end
            if !isleaf(tree, node)
                Recurse!(tree, node.other[1].inout[2]);
                Recurse!(tree, node.other[2].inout[2]);
            end
        end
        
        root = getroot(tree)
        Recurse!(tree, root) # do the recursive simulations
        tree
    end # menuramat!
    
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
        end
        collect(Iterators.flatten(res))
    end # predictTraitTree
    
    ##########################################################################33
    ###########################################################################3
    
    function putp!(tree, p1, key)
        for i in 1:length(tree.nodes) ## 'other' means branches
            tree.nodes[i].data[key] = p1;
        end
        tree;
    end # putp!
    
    ####################################################################33
    ######################################################################3
    function gen_cov_mat(mat, p, tspan, u0=zeros(size(mat)), dt = 0.001)
        function drift(du, u, p, t) ## drift function for the SDE
            du .= p.a * t .* p.A
        end # drift
        
        function diffusion(du, u, p, t) ## diffusion function for the SDE
            du .= p.b * t .* p.B 
        end # diffusion
        
        lowertri = LowerTriangular(mat);
        uppertri = - UpperTriangular(mat);
        skewsymm = lowertri + uppertri;
        pp = (A=skewsymm, B=skewsymm, a=p.a, b=p.b); ## skew symmetric matrices not necessarily the same.
        prob = SDEProblem(drift, diffusion, u0, tspan, p=pp); ## setup SDE problem
        sol = solve(prob, EM(), p=pp, dt=dt);
        Omega1 = last(sol.u); ## get the final matrix
        
        exp(Omega1) * mat * exp(-Omega1) ## reconstruct P_1
    end # gen_cov_mat
    ######################################################################3
    ##################################################################3n


    tr = Ultrametric(100);
    a1=1.0;
    b1=0.1;
    x0 =  [5.843333, 3.057333, 3.758000, 1.199333]
    tree = rand(tr); 
       time_tot = 1.0
    tspan = (0.0, time_tot)
    P0 =  [0.6856935  -0.0424340    1.2743154   0.5162707;
       -0.0424340   0.1899794   -0.3296564  -0.1216394;
       1.2743154  -0.3296564    3.1162779   1.2956094;
       0.5162707  -0.1216394    1.2956094   0.5810063];
alpha1 = (3.0, 3.0, 3.0, 3.0);
mu1 = (5.843333, 3.057333, 3.758000, 1.199333); ## Start at the trait means
sigma1 = (1.0, 1.0, 1.0, 1.0);

function simulate((alpha, mu, sigma))
    p1 = (alpha, mu, sigma, mat=P0, a=a1, b=b1)
    putp!(tree, p1, "parameters");
    menuramat!(tree);
    menura!(tree);
(tree, predictTraitTree(tree))
end

exampledat = simulate((alpha1, mu1, sigma1));

priordists = [Truncated(Normal(0,3), 0, Inf)]
priors = Factored(repeat(priordists, 12)...,)



function cost((alpha, mu, sigma))
    x=simulate((alpha, mu, sigma))[2]
    y=exampledat[2]
    euclidean(x, y)
end #cost

approx_density = ApproxKernelizedPosterior(priors,cost,0.005)
res = sample(approx_density, AIS(25), 10000, ntransitions=100, discard_initial = 250)   

save_object("ABCResults.jld2", res)

#  plot(tree,
#     size = (400, 800),
#     markersize = 20, 
#     series_annotations = text.(1:nnodes(tree), 15, :center, :center, :white))


# testnodes = getnodes(exampledat[1]);
#     plot(xlim = (0.0,1.0), ylim = (4.0, 8.0), zlim=(0.0, 5.0), legend=nothing,
#      reuse=false)
# for i in testnodes
#     u1= i.data["trace"]
#     uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
#     myt = i.data["timebase"]
#     plot!(myt, uu1[:, 1], uu1[:,2])
# end
# current()

# plot(xlim=(0.0,1.0), ylim= (2.0, 6.0), legend=nothing, reuse=false)
# for i in testnodes
#     u1= i.data["trace"]
#     uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
#     myt = i.data["timebase"]
#     plot!(myt, uu1[:, 3])
# end
# current()

plot(xlim=(2.0,6.0), ylim= (4.0, 8.0), legend=nothing, reuse=false)
for i in testnodes
    u1= i.data["trace"]
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    plot!(uu1[:,3], uu1[:, 1])
end
current()

plot(xlim=(2.0,6.0), ylim= (4.0, 8.0), zlim=(0.0, 5.0), legend=nothing,
     reuse=false)
for i in testnodes
    u1= i.data["trace"]
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["timebase"]
    plot!(uu1[:,3], uu1[:, 1], uu1[:,2])
end
current()

