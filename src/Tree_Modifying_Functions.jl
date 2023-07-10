using Phylo

######################################
# Functions to modify tree variables #
######################################

function putp!(tree, p1, key)
    for i in eachindex(tree.nodes)
        tree.nodes[i].data[key] = p1
    end
    tree;
end # putp!

function apply_trait(tree, trait, node::Int, key)
    tree.nodes[node].data[key] = trait
end


function apply_trait(tree, trait, node::Array{Int}, key)
    for n in node
        tree.nodes[n].data[key] = trait
    end
end

function apply_trait(tree, trait, node::LinkNode, key)
    node.data[key] = trait
end

function apply_trait_descend(tree, trait, node, key)
    tree.nodes[node].data[key] = trait
    for i in getdescendants(tree,tree.nodes[node])                                                                   
        apply_trait(tree, trait, i, key)                                                                                                                   
    end
end

function apply_prior(tree, prior, node)
    apply_trait(tree, prior, node, "Prior")
end

function apply_prior_descend(tree, prior, node)
    apply_trait_descend(tree, prior, node, "Prior")
end

function get_trait(tree, key, node)
    tree.nodes[node].data[key]
end