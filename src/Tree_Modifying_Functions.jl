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

function apply_trait(tree, trait, node, key)
    tree.nodes[node].data[key] = trait
end

function apply_prior(tree, prior, node)
    apply_trait(tree, prior, node, "Prior")
end