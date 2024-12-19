library(ape)
tee <- read.tree("pruned7.tre")

max(diag(vcv.phylo(tee)))

tee.size <- max(diag(vcv.phylo(tee)))

tee$edge.length <- tee$edge.length/tee.size

write.tree(tee, file = "prunedscaled.tre")
