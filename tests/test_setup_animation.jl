# Test file which sets up a basic instance for a run
using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

pyplot()

Random.seed!(1)

# number of parameters
n = 9

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)
time_tot = 1.0
tspan = (0.0, time_tot)

# G matrix
P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  ))) # Change to Wishart Distribution

# traits needed to evolve traits
alpha1 = repeat([1.0], n)
mu1 = repeat([0.0], n) ## randn(8); ## Start at the trait means
sigma1 = repeat([1.0], n)
trait_para = (alpha=alpha1, sigma=sigma1, mu = mu1)

# Variables needed for OU matrix model
mat_alpha = 0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)
mat_para = (alpha = mat_alpha, sigma = mat_sigma, mu = mat_mu)

data = menura_sim_exper_mat_OU(alpha1, sigma1, mu1, P0, mat_alpha, mat_sigma, mat_mu, tree1)

# Plotting 

# Change backend
gr()

# Set up 2d points
points = Vector{Vector{Any}}()

# Create vector of time and trait values
for (col_i, node) in enumerate(data[1].nodes) 
    nodepoints = [(time, [node.data["trace"][i][1], node.data["trace"][i][2]]) 
                    for (i, time)  in enumerate(node.data["timebase"])]
    push!(points, nodepoints)
end

s(x) = x[1][1]

sorted_points= sort!(points, by = s)

# create animation
anim = @animate for cut_off in 0.0:0.005:1.0
    xss = Vector()
    yss = Vector()
    for nodepoints in sorted_points 
        current_index = 1
        cut_off < nodepoints[current_index][1] && continue # check if first element too large
        for nodepoint in nodepoints # Get index up to required value
            nodepoint[1] < cut_off && (current_index += 1)
        end
        current_index -= 1
        xs = [nodepoints[i][2][1] for i in 1:current_index]
        ys = [nodepoints[i][2][2] for i in 1:current_index]
        
        push!(xss, xs)
        push!(yss, ys)
    end
    p = plot(xss, yss, xlim = (-2.5,2.5), ylim = (-2.5, 2.5), legend=nothing, reuse=false)
end

# create gif for 2d animation
gif(anim, "Evoling_traits.gif", fps = 20)

# set up points for 3d animation
points3d = Vector{Vector{Tuple{Float64, Vector{Float64}}}}()
for node in data[1].nodes 
    nodepoints = [(time, [node.data["trace"][i][1], node.data["trace"][i][2], 
                            node.data["trace"][i][3]]) 
                            for (i, time)  in enumerate(node.data["timebase"])]
    push!(points3d, nodepoints)
end

# sort points by starting time. This ensures the colours are consistent
s(x) = x[1][1]

sorted_points3d = sort!(points3d, by = s)

# create 3d animation
anim = @animate for cut_off in 0.0:0.005:1.0
    xss = Vector()
    yss = Vector()
    zss = Vector()
    for nodepoints in sorted_points3d 
        current_index = 1
        cut_off < nodepoints[current_index][1] && continue # check if first element too large
        for nodepoint in nodepoints # Get index up to required value
            nodepoint[1] < cut_off && (current_index += 1)
        end
        current_index -= 1
        xs = [nodepoints[i][2][1] for i in 1:current_index]
        ys = [nodepoints[i][2][2] for i in 1:current_index]
        zs = [nodepoints[i][2][3] for i in 1:current_index]
        push!(xss, xs)
        push!(yss, ys)
        push!(zss, zs)
    end
    p = plot(xss, yss, zss, xlim = (-2.5,2.5), ylim = (-2.5, 2.5), 
                zlim = (-2.5, 2.5),legend=nothing, reuse=false)
    plot!(p[1], camera = (45 + 30 * (sin(2*cut_off*π)), 40))
end

# create gif of 3d animation
gif(anim, "Evoling_traits_3d.gif", fps = 20)

# set up points for through time animation
points_through = Vector{Vector{Vector{Float64}}}()
for node in data[1].nodes 
    nodepoints = [[node.data["timebase"][i], node.data["trace"][i][1], 
                            node.data["trace"][i][2]] 
                            for (i, time)  in enumerate(node.data["timebase"])]
    push!(points_through, nodepoints)
end

# sort points by starting time. This ensures the colours are consistent
s2(x) = x[1]

sorted_points_through = sort!(points_through, by = s2)

# create through time animation
anim = @animate for cut_off in 0.0:0.005:1.0
    tss = Vector()
    xss = Vector()
    yss = Vector()
    for nodepoints in sorted_points_through 
        current_index = 1
        cut_off < nodepoints[current_index][1] && continue # check if first element too large
        for nodepoint in nodepoints # Get index up to required value
            nodepoint[1] < cut_off && (current_index += 1)
        end
        current_index -= 1
        xs = [nodepoints[i][2] for i in 1:current_index]
        ys = [nodepoints[i][3] for i in 1:current_index]
        ts = [nodepoints[i][1] for i in 1:current_index]
        push!(xss, xs)
        push!(yss, ys)
        push!(tss, ts)
    end
    p = plot(tss, xss, yss, xlim = (0,1), ylim = (-2.5, 2.5), 
                zlim = (-2.5, 2.5),legend=nothing, reuse=false)
end

# create gif of through time animation
gif(anim, "Evoling_traits_through.gif", fps = 20)