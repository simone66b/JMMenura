using DifferentialEquations, Phylo, Plots, Distributions, StatsPlots,
Distances, KissABC, LinearAlgebra, JLD2; pyplot();
## import Phylo.API;

## Iris data:

mat = [1.0000000  -0.1175698    0.8717538   0.8179411;
      -0.1175698   1.0000000   -0.4284401  -0.3661259;
       0.8717538  -0.4284401    1.0000000   0.9628654;
       0.8179411  -0.3661259    0.9628654   1.0000000]

## mat=Matrix(1.0I, 4, 4) 

time_tot = 1.0
tspan = (0.0, time_tot)
x0 = [5.843333, 3.057333, 3.758000, 1.199333] ## actual iris means
## x0 = [0.8, 0.3, 0.375, 0.199]
## tree = open(parsenewick, "tree.phy") Use if reading tree in from file
## We will create a random tree instead.

alpha_vec = [3.0, 3.0, 3.0, 3.0]
mu_vec = x0 ##  [0.5, 0.43, 0.5, 0.5]
sigma_vec = [1.0, 1.0, 1.0, 1.0]

## x0 = [0.5, 0.5, 0.5, 0.5]
p = [alpha_vec, mu_vec, sigma_vec]

function simulation(x0, mat, tspan, p)
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

function menura(tree, x0, mat, p, t0 = 0.0, dt=0.001) 
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
end # menura

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


    tr = Ultrametric(100);
    tree = rand(tr);
    plot(tree)
    tree2  = menura(tree, x0, mat, p)
    predictTraitTree(tree2);
    tree2
end # simulation

### plot(tree)
 test = simulation(x0, mat, tspan, p)


## euclidean(vals, resvalsflat)
###########################3 Plotting ###################################3
testbranches = getbranches(test);
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

