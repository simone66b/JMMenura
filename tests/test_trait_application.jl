using .JMMenura
using Phylo, Distributions, Random, Plots, LinearAlgebra

# Testing how traits are applied.

Random.seed!(1)

# Creating tree
tree = Ultrametric(6)
Random.seed!(1)
tree1 = rand(tree)

# viewing tree
plot_labelled(tree1)

# setting up Variables
mat_parameters = Dict(11 => (alpha = 1.0, mu = 1.0, sigma = 1.0), 8 => (alpha = 2.0, mu = 2.0, sigma = 2.0))
trait_parameters = Dict(11 => (alpha = 2.0, mu = 3.0, sigma = 1.0), 8 => (alpha = 6.0, mu = 0.0, sigma = 20.0))


mat_keys = sort(collect(keys(mat_parameters)), rev = true)
trait_keys = sort(collect(keys(trait_parameters)), rev = true)

for mat_key in mat_keys 
    apply_trait_descend(tree1, mat_parameters[mat_key], mat_key, "mat_para")
end

for trait_key in trait_keys 
    apply_trait_descend(tree1, trait_parameters[trait_key], trait_key, "trait_para")
end

# Check to see if worked

for node in getnodes(tree1)
    println(node)
    println("Trait parameters: $(node.data["trait_para"])")
    println("Matrix parameters: $(node.data["mat_para"])")
    println()
end