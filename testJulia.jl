println("hello world")
using Phylo, Random, Distributions

tree = Ultrametric(100)
using Plots, Random
plot(tree)
rand!(BrownianTrait(tree, "Trait"), tree)
plot(tree, line_z="Trait")
     
using DiffEqFlux, DifferentialEquations, Plots
function BM!(du,u,t)
  x = u
  du = dx = 0.0
end
u0 = 0.0
tspan = (0.0,1.0)

function noise!(du,u,t)
  x = u
  du = 0.3
end
#= p = [1.5,1.0,3.0,1.0,0.3,0.3] =#

prob = SDEProblem(BM!, noise!,u0,tspan)
sol = solve(prob) 
plot(sol)


#= ############################################################ =#

using StochasticDiffEq, StaticArrays
const μ = 0.5ones(4)
const σ = 0.1ones(4)
f(du,u,p,t) = du .= μ .* u
g(du,u,p,t) = du .= σ .* u
u0 = 0.1ones(4)
tspan = (0.0,1.0)
saveat = range(0.0,1.0,length=20)
prob = SDEProblem(f,g,u0,tspan)

@time for i in 1:100
  sol = solve(prob,SRIW1(),adaptive=false,dt=saveat[2])
end

using StochasticDiffEq, StaticArrays
mu = 0.0
sigma = 0.1
alpha = 0.2
f(du,u,p,t) = du .= alpha * (mu .- u) + sigma .* u
u0 = 0.1
tspan = (0.0,1.0)
saveat = range(0.0,1.0,length=1000)
prob = SDEProblem(f,u0,tspan)

@time for i in 1:100
  sol = solve(prob,SRIW1(),adaptive=false,dt=saveat[2])
end

#= ######################################################### =#
using DifferentialEquations, Plots

function drift(dx,x,p,t)
    alpha, mu, sigma, tau = p; 
    dx = alpha * (mu - x)
end


function stoch(dx,x,p,t);
  alpha, mu, sigma, tau = p; 
  dx = sqrt((2.0*sigma^2.0/tau))
end

x0 = 0.0
p = (0.3, 0.0, 1.5, 0.01);
tspan=(0.0, 50.0);
prob = SDEProblem(drift, stoch, x0, tspan, p)
sol = solve(prob); #= lots of errors here!! =#


#= ################################################################3 =#

using DifferentialEquations, Plots, DifferentialEquations.EnsembleAnalysis 

    
function diff(x,p, t)
    mu, sigma, alpha = p;
    sigma
end

function drift(x, p, t)
    mu, sigma, alpha = p;
    alpha * (mu - x)
end

dt = 0.005
x0 = 0.0
p = (0.1, 1.0, 1.5);
tspan = (0.0, 1.0);
prob = SDEProblem(drift, diff, x0, tspan,p)
@time begin
    sol =  solve(prob,EM(), dt=dt)
end
#= 0.006506  seconds =#
plot(sol)


ensembleprob = EnsembleProblem(prob)
sol = solve(ensembleprob,EnsembleThreads(),trajectories=1000)
summ = EnsembleSummary(sol, 0.0:0.01:1.0)
plot(summ,labels="Middle 95%")

summ = EnsembleSummary(sol,0:0.01:1;quantiles=[0.25,0.75])
plot!(summ,labels="Middle 50%",legend=true)


#= ###############################################################33 =#
using Random, DifferentialEquations, Plots

mutable struct Node
    name::Int64
    parent::Union{Int64, Missing}
    left::Union{Int64, String}
    right::Union{Int64, String}
    llength::Float64
    rlength::Float64
    lphenotype::Any
    rphenotype::Any
    rheight::Union{Float64, Missing}
    lheight::Union{Float64, Missing}
end

tree = ["z", 1.0, [["a", 1.0, ["b", 1.0, "c", 1.0], 1.0], 1.0, ["d", 1.0, "e", 1.0], 1.0], 1.0]
tree = [["a", 1.0, "b", 1.0], 1.0, ["c", 1.0, "d", 1.0], 1.0]
tree = ["a", 1.0, "b", 1.0]


function unfold(tree)
    treelist =  []
    index = 1
    function recurse(tree, index)
        if isa(tree[1], String) && isa(tree[3], String)
            theNode =  Node(index, missing, tree[1] , tree[3], tree[2], tree[4], [], [], missing, missing) 
            push!(treelist, theNode)
            index+=1
        elseif isa(tree[1], Vector) && isa(tree[3], String)
            theNode = Node(index, missing, index+=1, tree[3], tree[2], tree[4], [], [], missing, missing)
            push!(treelist, theNode)
            recurse(tree[1], index)
        elseif isa(tree[1], String) && isa(tree[3], Vector)
            theNode = Node(index, missing, tree[1], index+=1, tree[2], tree[4], [], [], missing, missing)
            push!(treelist, theNode)
            recurse(tree[3], index)
        elseif  isa(tree[1], Vector) && isa(tree[3], Vector) 
            theNode = Node(index, missing, index+=1, index+=1, tree[2], tree[4], [], [], missing, missing)
            push!(treelist, theNode)
            recurse(tree[1], index)
            recurse(tree[3], index)
        end
    end
    recurse(tree, index)
    treelist
end

treelist = unfold(tree)
##################################################################

using DifferentialEquations, Plots
mutable struct Node
    name::Int64
    parent::Union{Int64, Missing}
    left::Union{Int64, String}
    right::Union{Int64, String}
    llength::Float64
    rlength::Float64
    lphenotype::Any
    rphenotype::Any
    rheight::Union{Float64, Missing}
    lheight::Union{Float64, Missing}
end



Node1 = Node(1, missing, 2, 3, 4.5, 4.5, [], [], missing, missing)
Node2 = Node(2, 1, "a", 4, 2.0, 2.0, [], [], missing, missing)
Node3 = Node(3, 1, "b", "c", 2.0, 2.0, [], [], missing, missing)
Node4 = Node(4, 2, "d", "e", 2.0, 2.0, [], [], missing, missing)
tree = [Node1, Node2, Node3, Node4]

function mysim(tree)
    nodeInit = findall(x -> ismissing(x.parent) == 1, tree)[1] # find the root
    ## define the OU simulation
    function OU(x0, tspan)
        function drift(x, p, t)
            mu, sigma, alpha = p;
            alpha * (mu - x)
        end
        
        function diff(x,p, t)
            mu, sigma, alpha = p;
            sigma
        end

        dt = 0.005
        p = (1, 0.3, 1.5);
        prob = SDEProblem(drift, diff, x0, tspan,p)
        sol =  solve(prob, EM(), dt=dt)
    end

    
        function Heights!(tree, nodeInit)
    ## get the node heights so we know when to start the simulations
            node = tree[nodeInit]
            if ismissing(node.parent)
                node.lheight = node.llength
               node.rheight = node.rlength
            else node.lheight=node.llength + tree[node.parent].lheight
                node.rheight = node.rlength + tree[node.parent].rheight
            end
            if isa(node.left, Int64) Heights!(tree, node.left)
            end
            if isa(node.right, Int64) Heights!(tree, node.right)
            end
        end
            
    function Recurse!(tree, nodeInit)
        ## recurse along the tree, simulating the OU process for each branch and storing it in lphenotype and rphenotype
        node = tree[nodeInit]

        if ismissing(node.parent) ## the root node, to get started
            node.lphenotype = OU(0.0, (0.0, node.lheight))
            node.rphenotype = OU(0.0, (0.0, node.rheight))
            
        elseif tree[node.parent].left == node.name
            node.lphenotype = OU(tree[node.parent].lphenotype[end],
                                        (tree[node.parent].lheight, node.lheight))
                    
            node.rphenotype = OU(tree[node.parent].lphenotype[end],
                                        (tree[node.parent].rheight, node.rheight))
        elseif tree[node.parent].right == node.name
            node.lphenotype = OU(tree[node.parent].rphenotype[end],
                                        (tree[node.parent].lheight, node.lheight))
            node.rphenotype = OU(tree[node.parent].rphenotype[end],
                                        (tree[node.parent].rheight, node.rheight))
        end
        if isa(node.left, Int64) Recurse!(tree, node.left) 
        end
        if isa(node.right, Int64) Recurse!(tree, node.right)
        end
    end
    Heights!(tree, nodeInit) # calculate the heights
    Recurse!(tree, nodeInit) # do the recursive simulations
    tree
end

test = mysim(tree)

plot(xlim=[0.0, 9.0], ylim= [0.0, 1.5], legend=nothing)
for i in 1:length(test)
    plot!(test[i].lphenotype.t, test[i].lphenotype.u, legend= nothing)
    plot!(test[i].rphenotype.t, test[i].rphenotype.u, legend = nothing)
end
current()
savefig("test.png")


#####################################################################33

using Phylo, Plots

tree = parsenewick("((a:1.0, (b:1.0, c:1.0):1.0):1.0, d:1.0);")
theRoot = getroot(tree)

theRoot.other[1].data["1"] = rand(10)
