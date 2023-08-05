using Phylo
using .JMMenura
using Distributions

tree = Ultrametric(9)
tree1 = rand(tree)

plot_labelled(tree1)

# apply_prior_descend(tree1, Exponential(1), 16)
