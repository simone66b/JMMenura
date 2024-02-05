using PosDefManifold, LinearAlgebra, DifferentialEquations, StatsBase, Plots
using CSV, Tables, DataFrames, XLSX, Statistics, Pipe, Phylo
pyplot();

cd("/home/simoneb/Desktop/JMMenura") ## may need to change this to suit

#############################################################################
function drift(du, u, p, t)
    du .= p.alpha .* logMap(Fisher, Hermitian(u), Hermitian(p.mu))# Mean reversion
end  ## drift function

function diffusion(du, u, p, t)
    du .= p.sigma  ## scaled BM
end ## diffusion function

##function OUmatrix(drift, diffusion, G0, mu, alpha, sigma, tspan, dt)
temp = Matrix{Float64}[]
push!(temp, G0)
pp = (mu = mu, G0=G0, alpha =  alpha, sigma = sigma) ## tuple of parameters

prob = SDEProblem(drift, diffusion, G0, tspan, p = pp) ## set up the sde problem
sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
G1 = expMap(Fisher, Hermitian(sol.u[2]), Hermitian(sol.u[1])) # back transfer trace onto manifold
push!(temp, G1)
  ## can now plot, print this, etc.
end # OUmatrix function




###################################################################################################3

nms = ["cris.txt", "ever.txt", "grah.txt", "line.txt", "pulc.txt", "sagr.txt", "smar.txt"];
matmatrix = CSV.File(nms, header = 0) |> Tables.matrix
Garray = reshape(matmatrix', 8, 8, 7) ## 8 traits, 7 species
Gvec = [Garray[1:8, 1:8, i] for i in 1:7] ## Vector of Species Gmatrices
DFmats = DataFrame(Gmatrix = Gvec)

xf = DataFrame(XLSX.readtable("Adult measurements for divergence.xlsx", "Pmatrix Measurements with outli"))

colNames = append!([2], 4:11) # columns to keep
sppNames = ["A. cristatellus", "A. evermanni", "A. grahami", "A. lineatopus",
 "A. pulchellus", "A. sagrei", "A. smaragdinus"]
AnolesData = @pipe xf[:, colNames] |> 
filter(row -> row.Species in sppNames)|>
groupby(_, :Species) |>
combine(_, 2:9 .=> (x -> mean(skipmissing(x)))) |> ## species means
transform(_, 2:9 .=> (x -> log.(x))) |> ## take logs
select(_, Not(2:9)) |> 
rename(_, names(xf)[colNames])

FullData = hcat(AnolesData, DFmats)

tree = open(parsenewick, Phylo.path("/home/simoneb/Desktop/JMMenura/pruned7.tre"))
## display(plot(tree, linewidth=5, tipfont=20))

G0Vec = [0.285, 0.118, 0.131, 0.053, 0.212, 0.172, 0.188, 0.210, 0.277,
            0.230, 0.103, 0.074, 0.073, 0.045, 0.023, 0.800, 0.317, 0.054,
            0.121, 0.152, 0.147, 0.654, 0.028, 0.047, 0.124, 0.148, 0.858,
            0.749, 0.581, 0.626, 0.906, 0.553, 0.644, 0.708, 0.731, 0.872]

            dt = 0.001
G0 =zeros(8, 8)
idx = tril!(trues(size(G0)))
G0[idx] .= G0Vec
G0 .= Hermitian(G0, :L)
## u0 = G0 #### FullData[7, :Gmatrix]; ## randPosDefMat(8)
mu = G0
alpha = fill(2.0, size(G0, 1), size(G0, 1));
sigma = fill(1 / sqrt(2), size(G0, 1), size(G0, 1));
sigma[diagind(sigma)] .= 1;
tspan = (0.0, 0.001)

@time test = OUmatrix(drift, diffusion, G0, mu, alpha, sigma, tspan, dt);
#= # Extract the elements from the matrices
test2 = convert(Vector{Matrix{Float64}}, test);
## test2 .= cov2cor.(test2)  ## uncomment to see correlation matrix plots
vals = map(x -> x[tril(trues(size(x)))], test2); 
plot(reuse = false) ## extract lower triangular from test
for i in 1:Integer(dim(test[1], 1) * (dim(test[1], 1) - 1) / 2 + dim(test[1], 1))
plot!(map(x -> x[i], vals), legend = false) =#
end
## current()
