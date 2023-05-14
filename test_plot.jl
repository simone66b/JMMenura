using Phylo
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


tree = Ultrametric(39)
tree1 = rand(tree)
sort!(tree1)
p = plot(tree1)
show(p)
a = round.(nodeheights(tree1), sigdigits = 10)
b = length_order!(a, tree1)

plot(tree1, linewidth = 5,
           markersize = 25,
           series_annotations = text.(b, 5, :center, :center, :white))
