println("hello world")
using Phylo, Random, Distributions

tree = rand(Ultrametric(2048, 10.0))
using Plots, Random
plot(tree)
@time rand!(BrownianTrait(tree, "Trait"), tree)
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

using DifferentialEquations, Phylo, Plots; pyplot()

function mysim(tree, x0=0.0, t0 = 0.0)
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

tr = Ultrametric(1024, 10.0);
tree = rand(tr);
plot(tree)
savefig("tree.png")

test = mysim(tree);
testbranches=getbranches(test);

plot(xlim=[0.0, 10.0], ylim= [-0.1, 2.0], legend=nothing)
for i in testbranches
    plot!(i.data["1"].t, i.data["1"].u, legend= nothing)
end
current()
savefig("trace.png")


#####################################################################33

using Phylo, Plots

tree = parsenewick("((a:1.0, (b:1.0, c:1.0):1.0):1.0, d:1.0);")
theRoot = getroot(tree)

theRoot.other[1].data["1"] = rand(10)

#######################################################################
## import Pkg; Pkg.add("PyPlot")

using DifferentialEquations, Phylo, Plots, Distributions, StatsPlots; pyplot()

alpha_vec = [3.0, 3.0, 4.0, 2.0]

mu_vec = [5.843333, 3.057333, 3.758000, 1.199333]
sigma_vec = [1.0, 1.0, 1.0,  1.0]

p = [alpha_vec, mu_vec, sigma_vec]
function f(du, u, p, t)
    alpha , mu, sigma = p
    du .= alpha .* (mu - u)
end # f

function g(du, u, p, t)
    alpha, mu, sigma = p
    du .= sigma
end # g

# rho = 0.0
# mat = [1.0 rho rho rho rho;
#        rho 1.0 rho rho rho;
#        rho rho 1.0 rho rho;
#        rho rho rho 1.0 rho;
#        rho rho rho rho 1.0]

mat = [1.0000000  -0.1175698    0.8717538   0.8179411;
      -0.1175698   1.0000000   -0.4284401  -0.3661259;
       0.8717538  -0.4284401    1.0000000   0.9628654;
       0.8179411  -0.3661259    0.9628654   1.0000000]

time_tot = 100.0
tspan = (0.0, time_tot)

OU_noise = CorrelatedWienerProcess(mat, tspan[1], zeros(4), zeros(4)) 

u0 = mu_vec

prob = SDEProblem(f,g, u0,tspan, p=p, noise=OU_noise)
dt = 1e-3

nums = Int(time_tot/dt) + 1
sol = solve(prob, dt=dt, p=p, adaptive=false);
plot(sol)

myt = sol.t
u1 = sol.u

uu1 = transpose(reshape(collect(Iterators.flatten(u1)), 4, nums))

plot(myt, uu1[1:nums,3], uu1[1:nums,4])
plot(uu1[1:nums,1], uu1[1:nums,4])

histogram(uu1[1:nums, 1])
plot(qqplot(Normal, uu1[1:nums, 1]))
