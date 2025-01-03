library(ape)
tee <- rtree(9)

tee=chronos(tee, lambda=0)  

write.tree(tee, file = "sim.tre")
