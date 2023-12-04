
##########################################
# Functions which Handle Simulation Loop #
##########################################

"""
recurse_menura!(tree, node, trait_evol, matrix_evol, each::Bool)
Recursively iterates over the tree and evolves each node.

Takes in tree, node to evol and functions specifing how to evolve.
"""
function recurse_menura!(tree, node, trait_evol::Function, matrix_evol::Function, each::Bool)
    if ismissing(node.inbound) ## the root node, to get started
        # Functionality now to be done by a higher up function.
        # node.data["mat_trace"] = node.data["parameters"].mat ## starting matrix

        # node.data["trait_trace"]= [x0]
        # node.data["timebase"] = [t0]
    else    
        ancestor = getancestors(tree, node)[1]

        if haskey(ancestor.data, "known_G_mat")
            evol_matrix = ancestor.data["known_G_mat"] # should this be ancestor? I think so
         else 
             evol_matrix = ancestor.data["mat_trace"][end]
         end

        node.data["mat_trace"] =
            matrix_evol(evol_matrix,ancestor.data["mat_para"], 
                        (getheight(tree, ancestor),getheight(tree, node)), each)

        evol = trait_evol(ancestor.data["trait_trace"][end], node.data["mat_trace"], 
                            ancestor.data["trait_para"],
                         (getheight(tree, ancestor), getheight(tree, node)), each)
        node.data["trait_trace"] = evol.u
        node.data["timebase"] = evol.t
        
    end
    if !isleaf(tree, node)
        recurse_menura!(tree, node.other[1].inout[2], trait_evol, matrix_evol, each)
        recurse_menura!(tree, node.other[2].inout[2], trait_evol, matrix_evol, each)
    end
end # Recurse! 


#################################
# Simulation Handling Functions #
#################################

# The function which runs the loop. Takes in a tree and the two evolution functions.
# Should throw an error if tree not set up properly
function menura!(tree, trait_evol::Function, matrix_evol::Function, t0::Float64, trait0::Vector{Float64}
                        , mat0::Matrix{Float64}, each::Bool)
    # set up starting time and values                     
    root = getroot(tree)

    root.data["trait_trace"]= [trait0]
    root.data["mat_trace"]= [mat0]
    root.data["timebase"] = [t0]

    # Getting cranky error function
    
    # perform the recursion
    recurse_menura!(tree, root, trait_evol, matrix_evol, each)

    # return the tree along with the final trait data points
    # Might have to reconsider this function
    (tree, reduce(hcat, [tip.data["trait_trace"][end] for tip in getleaves(tree)]))
end

# this function should now be a more fancy menura!
# Takes traits in dictionary and applies them in a descending pattern
function menura_parameter_descend!(mat_parameters, trait_parameters, tree, 
                                trait_evol::Function, matrix_evol::Function, t0::Float64, trait0::Vector{Float64}
                                , mat0::Matrix{Float64}, each::Bool)
# apply traits. Note phylo numbers the nodes in descending order by time
mat_keys = sort(collect(keys(mat_parameters)), rev = true)
trait_keys = sort(collect(keys(trait_parameters)), rev = true)

for mat_key in mat_keys 
    apply_trait_descend(tree, mat_parameters[mat_key], mat_key, "mat_para")
end

for trait_key in trait_keys 
    apply_trait_descend(tree, trait_parameters[trait_key], trait_key, "trait_para")
end

# run menura!
menura!(tree, trait_evol, matrix_evol, t0, trait0, mat0, each)

end
