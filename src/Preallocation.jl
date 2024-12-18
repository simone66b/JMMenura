

function preallocate_tree!(tree, dt, n)
    root = getroot(tree)
    children = getchildren(tree, root)

    recurse_preallocate_tree!(tree, children[1], dt, n)
    recurse_preallocate_tree!(tree, children[2], dt, n)
end

function recurse_preallocate_tree!(tree, node, dt, n)
    ancestor = getancestors(tree, node)[1]

    mat_num = floor((getheight(tree, node) - getheight(tree, ancestor))/dt) + 2

    node.data["trait_trace"] = [zeros(n) for _ in 1:mat_num]

    node.data["mat_trace"] = [zeros(n,n) for _ in 1:mat_num]

    node.data["timebase"] = [0.0 for _ in 1:mat_num]

    if !isleaf(tree, node)
        children = getchildren(tree, node)
        recurse_preallocate_tree!(tree, children[1], dt, n)
        recurse_preallocate_tree!(tree, children[2], dt, n)
    end
end