using DifferentialEquations, Phylo, Plots, Distributions, Distances, KissABC, JLD2, LinearAlgebra, ExponentialUtilities; pyplot();

function simulation(p1, tree)
    function diffusion(x0, tspan, p, dt=0.001)
        function drift(du, u, p, t)
            alpha = p[1];
            mu = p[2];
            ## du .= alpha .* u would be BM with drift
            du .= alpha .* (mu .- u) ## OU-like
        end # drift

        function diffusion(du, u, p, t)
            sigma = p[3]
            du .= sigma ## would be OU
            ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
            ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
        end # diff
        cov1=p[4]
        cor1 = cor(p[4])
        noise = CorrelatedWienerProcess(cor1, tspan[1],
                                        zeros(dim(cor1)),
                                        zeros(dim(cor1)))
      
    prob = SDEProblem(drift, diffusion, x0, tspan, p=p, noise=noise)       
    solve(prob, dt=dt, p=p, adaptive=false)
end # diffusion

function menura!(tree, x0, t0 = 0.0, dt=0.001) 
    function Recurse!(tree, node)
        if ismissing(node.inbound) ## the root node, to get started
            for i in 1:length(node.other) ## 'other' means branches
                node.other[i].data["1"]  = ## "1" is a placeholder for the Dict
                diffusion(x0, (t0, getheight(tree, node.other[i].inout[2])),
                          node.data["2"][i]);
            end
        else
            for i in 1:length(node.other) ## NB 'other' are branches here!
                node.other[i].data["1"] =
                    diffusion(node.inbound.data["1"].u[end],
                       (getheight(tree, node), getheight(tree, node) +
                        node.other[i].length),  node.data["2"][i]);
            end # for
        end
        for i in 1:length(node.other)
            if !isleaf(tree, node.other[i].inout[2])
                Recurse!(tree, node.other[i].inout[2]);
            end
        end
    end# Recurse!

    nodeInit = getroot(tree)
    Recurse!(tree, nodeInit) # do the recursive simulations
    tree
end # menura!

function predictTraitTree(tree)

    #### get the last multivariate trait value in a branch
    testtips = collect(nodefilter(isleaf, tree));
    res = Array{Vector{Float64}}(undef, length(testtips)); 
    tipnames = Array{String}(undef, length(testtips));
    tiptimes = Vector{Float64}();
   ##  finaltraitvals = Dict();
    ## tiptimesdict = Dict();
    
    for i in 1:length(testtips) ## could maybe use heightstoroot() for this computation
        res[i] = testtips[i].inbound.data["1"].u[end];
        tipnames[i] = testtips[i].name;
        push!(tiptimes, testtips[i].inbound.data["1"].t[end]);
    end

    collect(Iterators.flatten(res));
end # predictTraitTree

    function putp!(tree, p1, p2=p1)
        for i in 1:length(tree.nodes) ## 'other' means branches
        tree.nodes[i].data["2"] = (p1, p2);
        end
        tree;
    end # putp!

    ####################################################################33
    ######################################################################3
    function gen_cov_mat(p, tspan, dt = 0.001)
         function drift(du, u, p, t) ## drift function for the SDE
            du .= 1.0 * t .* p.A ## a = 1.0 for example
        end # drift

        function diffusion(du, u, p, t) ## diffusion function for the SDE
            du .= 1.0 * t .* p.B ## diffusion function b= 1.0 for example
        end # diffusion
        
        lowertri = LowerTriangular(p)
        uppertri = - UpperTriangular(p)
        skewsymm = lowertri + uppertri
        u0 = zeros(size(p))
        p=(A=skewsymm, B=skewsymm) ## skew symmetric matrices not necessarily the same.
        prob = SDEProblem(drift, diffusion, u0, tspan, p=p); ## setup SDE proble
        sol = solve(prob, ISSEM(theta=1/2, symplectic=true), p=p, dt=dt);
        Omega1 = last(sol.u); ## get the final matrix
        Omega2 = -last(sol.u); ## get another copy NB negative
        First = exponential!(Omega1); ## matrix exponential
        Second = exponential!(Omega2); ## matrix exponential
        result = First * p * Second; ## reconstruct P_1
        result
    end # gen_cov_mat
    
        ######################################################################3
        ##################################################################3
    tree = putp!(tree, p)
    tree2  = menura!(tree, x0)
   (tree2, predictTraitTree(tree2))
    ## tree2
end # simulation

    time_tot = 1.0
    tspan = (0.0, time_tot)
    x0 = [5.843333, 3.057333, 3.758000, 1.199333] ## starting values
P0 =  [0.6856935  -0.0424340    1.2743154   0.5162707;
       -0.0424340   0.1899794   -0.3296564  -0.1216394;
       1.2743154  -0.3296564    3.1162779   1.2956094;
       0.5162707  -0.1216394    1.2956094   0.5810063];


tr = Ultrametric(100)
tree = rand(tr)
 ## Q matrix
a1=1.0;
b1=2.0;

alpha1 = (3.0, 3.0, 3.0, 3.0)
mu1 = (5.843333, 3.057333, 3.758000, 1.199333); ## Start at the trait means
sigma1 = (1.0, 1.0, 1.0, 1.0);
p = (alpha1, mu1, sigma1, P0, a1, b1)
npar = 3 ## alpha mu, sigma of the SDE
ndims = length(alpha1)
exampledat = simulation(p, tree)







################################ ABC #################################################3
# priors = Factored(
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf),
#     Truncated(Normal(0, 3), 0, Inf))

priordists = [Truncated(Normal(0,3), 0, Inf)]
priors = Factored(repeat(priordists, 12)...,)



function cost((alpha1, mu1, sigma1, cov))
    x=simulation((alpha1, mu1, sigma1, cov), tree)[2]
    y=exampledat[2]
    euclidean(x, y)
end #cost

## approx_density = ApproxKernelizedPosterior(priors,cost,0.005)
## res = sample(approx_density, AIS(25), 10000, ntransitions=100, discard_initial = 250)   

## save_object("ABCResults.jld2", res)

## res = sample(approx_density,AIS(25), MCMCThreads(), 5000, 4, ntransitions=10, discard_initial = 250)
### ressmc = smc(priors, cost, nparticles=500, epstol=0.01)

 plot(tree,
    size = (400, 800),
    markersize = 20, 
    series_annotations = text.(1:nnodes(tree), 15, :center, :center, :white))


testbranches = getbranches(exampledat[1]);
    plot(xlim = (0.0,1.0), ylim = (4.0, 8.0), zlim=(0.0, 5.0), legend=nothing,
     reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(myt, uu1[:, 1], uu1[:,2])
end
current()

plot(xlim=(0.0,1.0), ylim= (2.0, 6.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(myt, uu1[:, 3])
end
current()

plot(xlim=(2.0,6.0), ylim= (4.0, 8.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    plot!(uu1[:,3], uu1[:, 1])
end
current()

plot(xlim=(2.0,6.0), ylim= (4.0, 8.0), zlim=(0.0, 5.0), legend=nothing,
     reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(uu1[:,3], uu1[:, 1], uu1[:,2])
end
current()

