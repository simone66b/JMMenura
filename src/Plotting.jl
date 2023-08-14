using Phylo, Plots
# Contains functions which allow for more easy plotting

function length_order_recurs!(array, min_height, label, iter)
    if iter > 0
        current_height, current_index = 1, nothing
        for (index,height) in enumerate(array)
            if height <= current_height && height >= min_height
                current_height, current_index = height, index
            end
        end
        if current_height == 1
            array[current_index] = 0 
            length_order_recurs!(array, current_height, label -1, iter -1 )
        else
            array[current_index] = label
            length_order_recurs!(array, current_height, label -1, iter -1 )
        end
    end
end

function length_order!(array, tree)
    length = nnodes(tree)
    length_order_recurs!(array, 0, length, length)
    a = zeros(length)
    j = 1
    for i in array
        if i != 0
            a[j] = i
            j += 1
        end
    end
    return a    
end

"""
    plot_labelled(tree)

Plots a phylogenetic tree numbering the nodes according to the indexing system used.

Useful for indexing into a tree.
"""
function plot_labelled(tree)
    a = round.(nodeheights(tree), sigdigits = 10)
    b = length_order!(a, tree)

    Plots.plot(tree, linewidth = 5,
            markersize = 25,
            series_annotations = Plots.text.(b, 5, :center, :center, :white))

end

function plot_data(tree, trait1, trait2; time = true)
    pyplot()
    p = Plots.plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
    if time 
        for node in getnodes(tree)
            u1= node.data["trace"]
            uui1 = [point[trait1] for point in u1]
            uui2 = [point[trait2] for point in u1]
            myt = node.data["timebase"]
            plot!(myt, uui1, uui2)
        end; # for
        return p
    else 
        for node in getnodes(tree)
            u1= node.data["trace"]
            uui1 = [point[trait1] for point in u1]
            uui2 = [point[trait2] for point in u1]
            plot!(uui1, uui2)
        end; # for
        return p
    end
end

function plot_data(tree, trait1; time = true)
    pyplot()
    p = Plots.plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
    for node in getnodes(tree)
        u1= node.data["trace"]
        uui1 = [point[trait1] for point in u1]
        myt = node.data["timebase"]
        plot!(myt, uui1)
    end; # for
    return p
end

function plot_data(tree, trait1, trait2, trait3; time = true)
    pyplot()
    p = Plots.plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
    for node in getnodes(tree)
        u1= node.data["trace"]
        uui1 = [point[trait1] for point in u1]
        uui2 = [point[trait2] for point in u1]
        uui3 = [point[trait3] for point in u1]
        plot!(uui1, uui2, uui3)
    end; # for
    return p
end

function plot_g_mat_evol(tree, leaf_num, name; fps = 10)
    nodes = [reverse(getancestors(tree, tree.nodes[leaf_num]))..., tree.nodes[leaf_num]]
    cov_mats = [(convert.(Matrix{Float64},node.data["matrixes"])) for node in nodes] 
    #THIS WILL NEED TO BE CHANGED

    lower = minimum([minimum(minimum.(cov_mat)) for cov_mat in cov_mats])
    upper = maximum([maximum(maximum.(cov_mat)) for cov_mat in cov_mats])

    anim = @animate for (i, cov_mat) in enumerate(Iterators.flatten(cov_mats))
        heatmap(cov_mat, title = "$i", yflip = true, clims = (lower, upper))
    end every 10

    gif(anim, name, fps = fps)
end

# animations for traits through time