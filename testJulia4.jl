using DifferentialEquations, Phylo, Plots, Distributions, Distances, KissABC; pyplot();

function simulation(ams)
    function diffusion(x0, tspan, p, dt=0.001)
        function drift(du, u, p, t)
            alpha, mu, sigma = p
            ## du .= alpha .* u would be BM with drift
            du .= alpha .* (mu .- u) ## OU-like
        end # drift

        function diff(du, u, p, t)
            alpha, mu, sigma = p
            du .= sigma ## would be OU
            ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
            ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
        end # diff

     noise = CorrelatedWienerProcess(mat, tspan[1],
                                    zeros(dim(mat)),
                                    zeros(dim(mat))) 
    prob = SDEProblem(drift, diff, x0, tspan, p=p, noise=noise)       
    solve(prob, dt=dt, p=p, adaptive=false)
end # diffusion

function menura!(tree, x0, mat, p, t0 = 0.0, dt=0.001) 
    function Recurse!(tree, node)
        if ismissing(node.inbound) ## the root node, to get started
            for i in 1:length(node.other) ## 'other' means branches
                node.other[i].data["1"]  = ## "1" is a placeholder for the Dict
               diffusion(x0, (t0, getheight(tree, node.other[i].inout[2])), p);
            end
        else
            for i in 1:length(node.other) ## NB 'other' are branches here!
                node.other[i].data["1"] =
                    diffusion(node.inbound.data["1"].u[end],
                       (getheight(tree, node), getheight(tree, node) +
                        node.other[i].length),p);
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
    
    ## map((i, j) -> finaltraitvals[i] = j, tipnames, res);
    ## map((i, j) -> tiptimesdict[i] = j, tipnames, tiptimes);
    collect(Iterators.flatten(res))
end # predictTraitTree

   tree2  = menura!(tree, x0, mat, p)
   predictTraitTree(tree2);
end # simulation

## tree = open(parsenewick, "exampletree.phy")
   ##  tr = Ultrametric(100);
   ##  tree = rand(tr);
    ## plot(tree)
    mat = [1.0000000  -0.1175698    0.8717538   0.8179411;
           -0.1175698   1.0000000   -0.4284401  -0.3661259;
           0.8717538  -0.4284401    1.0000000   0.9628654;
           0.8179411  -0.3661259    0.9628654   1.0000000]
    time_tot = 1.0
    tspan = (0.0, time_tot)
    x0 = [5.843333, 3.057333, 3.758000, 1.199333] ##


tr = Ultrametric(200)
tree = rand(tr)

alpha1 = (3.0, 3.0, 3.0, 3.0)
mu1 = (5.843333, 3.057333, 3.758000, 1.199333); ## Start at the trait means
sigma1 = (1.0, 1.0, 1.0, 1.0);
p = (alpha1, mu1, sigma1)

tst = simulation(p)

npar = 3 ## alpha mu, sigma of the SDE
ndims = length(alpha1)

priors = Factored(
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
                  Truncated(Normal(0, 3), 0, Inf),
Truncated(Normal(0, 3), 0, Inf),
Truncated(Normal(0, 3), 0, Inf))

exampledat = simulation(p)

function cost((alpha1, mu1, sigma1))
    x=simulation((alpha1, mu1, sigma1))
    y=exampledat
    euclidean(x, y)
end #cost

approx_density = ApproxKernelizedPosterior(priors,cost,0.05)
res = sample(approx_density,AIS(25), MCMCThreads(), 5000, 4, ntransitions=10, discard_initial = 250)   


### ressmc = smc(priors, cost, nparticles=500, epstol=0.01)
