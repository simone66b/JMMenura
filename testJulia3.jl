using DifferentialEquations, Phylo, Plots, Distributions, StatsPlots, KissABC; pyplot();

## Iris data:

mat = [1.0000000  -0.1175698    0.8717538   0.8179411;
      -0.1175698   1.0000000   -0.4284401  -0.3661259;
       0.8717538  -0.4284401    1.0000000   0.9628654;
       0.8179411  -0.3661259    0.9628654   1.0000000]

time_tot = 1.0
tspan = (0.0, time_tot)
x0 = [5.843333, 3.057333, 3.758000, 1.199333] ## actual iris means

function mysim(tree, x0, mat, t0 = 0.0)
    ## define the OU simulation
    function OU(x0, tspan)
        function drift(du, u, p, t)
            alpha, mu, sigma = p
            ## du .= alpha .* u
            du .= alpha .* (mu - u)
        end # drift

        function diff(du, u, p, t)
            alpha, mu, sigma = p
            ## du .= sigma
            du .= sqrt.(u) .* sigma
        end # diff

        dt = 0.001
        alpha_vec = [3.0, 1.0, 2.0, 3.0]
        mu_vec = [5.0, 3.0, 2.0, 3.0]
        sigma_vec = [1.0, 1.0, 2.0, 3.0]

        p = [alpha_vec, mu_vec, sigma_vec]
        OU_noise = CorrelatedWienerProcess(mat, tspan[1],
                                           zeros(length(alpha_vec)), zeros(length(alpha_vec))) 

        prob = SDEProblem(drift, diff, x0,tspan, p=p, noise=OU_noise)       
        sol = solve(prob, dt=dt, p=p, adaptive=false)   
    end #OU
    
    function Recurse!(tree, node)
        if ismissing(node.inbound) ## the root node, to get started
            for i in 1:length(node.other) ## 'other' means branches
                node.other[i].data["1"]  = ## "1" is a placeholder for the Dict
                OU(x0, (t0, getheight(tree, node.other[i].inout[2])));
            end
        else
            for i in 1:length(node.other) ## NB 'other' are branches here!
                node.other[i].data["1"] =
                    OU(node.inbound.data["1"].u[end],
                       (getheight(tree, node), getheight(tree, node) +
                        node.other[i].length))
            end
        end
        for i in 1:length(node.other)
            if !isleaf(tree, node.other[i].inout[2])
                Recurse!(tree, node.other[i].inout[2])
            end
        end
    end# Recurse!
    nodeInit = getroot(tree)
    Recurse!(tree, nodeInit) # do the recursive simulations
    tree
end # mysim

## tree = open(parsenewick, "tree.phy") Use if reading tree in from file
## We will create a random tree instead.

tr = Ultrametric(100, 1.0);
tree = rand(tr);
display(plot(tree))

test = mysim(tree, x0, mat);
testbranches=getbranches(test);

plot(xlim=(0.0,1.0), ylim= (0.0, 10.0), zlim=(0.0, 10.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(myt, uu1[:, 1], uu1[:,2])
end
current()

plot(xlim=(0.0,1.0), ylim= (0.0, 10.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(myt, uu1[:, 3])
end
current()

plot(xlim=(0.0,10.0), ylim= (0.0, 10.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    plot!(uu1[:,3], uu1[:, 1])
end
current()

plot(xlim=(0.0,10.0), ylim= (0.0, 10.0), zlim=(0, 10.0), legend=nothing, reuse=false)
for i in testbranches
    u1= i.data["1"].u
    uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, length(u1)))
    myt = i.data["1"].t
    plot!(uu1[:,3], uu1[:, 1], uu1[:,2])
end
current()

## savefig("trace.png")
