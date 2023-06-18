using Phylo, Plots
# Contains functions which allow for more easy plotting
pyplot()

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


# testnodes = getnodes(exampledat[1])
# p2 = plot(xlim = (0.0,1.0), ylim = (-2.0, 2.0), zlim=(-2.0, 2.0), legend=nothing, reuse=false)
# for i in testnodes
#     u1= i.data["trace"]
#     uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)))
#     myt = i.data["timebase"]
#     plot!(myt, uu1[:, 1], uu1[:,3])
# end; # for
# current()
# display(p2)

# p3 = plot(xlim=(0.0,1.0), ylim= (-2.0, 2.0), legend=nothing, reuse=false)
# for i in testnodes
#     u1= i.data["trace"];
#     uu1 = transpose(reshape(collect(Iterators.flatten(u1)), n, length(u1)));
#     myt = i.data["timebase"];
#     plot!(myt, uu1[:, 3]);
# end; # for
# current()
# display(p3)