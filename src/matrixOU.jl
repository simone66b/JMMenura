using Pkg
using PosDefManifold, LinearAlgebra, DifferentialEquations, StatsBase, Plots
using CSV, Tables, DataFrames, XLSX, Statistics, Pipe, Phylo
pyplot();
Pkg.activate(".")
pwd()

# cd("/home/simoneb/Desktop/JMMenura") ## may need to change this to suit

#############################################################################
function drift(du, u, p, t)
    du .= p.alpha .* (p.mu - u) ## Mean reversion
end  ## drift function

function diffusion(du, u, p, t)
    du .= p.sigma ## scaled BM
end

function OUmatrix(uu0, mu1, alpha, sigma, tspan, dt)

# function vec2Hermitian(vec, size)
#     Hermitian([ i <= j ? vec[Integer(j * (j - 1)/ 2 + i)] : 0 for i = 1:size, j = 1:size ])
#    end

# alphavec = alpha[tril(trues(size(alpha)))]
# sigmavec = sigma[tril(trues(size(sigma)))]

# function logMapVec(x)
# @pipe x |> logMap(Fisher,Hermitian(_), HermId) |>
# _[tril(trues(size(_)))]
# end

HermId = Hermitian(1.0I(size(uu0, 1))) ## Identity matrix

## map mean matrix onto tangent space
uu1 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(uu0), HermId))
mu2 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(mu1), HermId))

pp = (mu = mu2, alpha =  alpha, sigma = sigma) ## tuple of parameters

prob = SDEProblem(drift, diffusion, uu1, tspan, p = pp) ## set up the sde problem
sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama

## herms = vec2Hermitian.(sol.u, size(uu0,1))
## herms = convert(Vector{Hermitian{Float64, Matrix{Float64}}}, herms)
temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
temp  ## can now plot, print this, etc.
end # OUmatrix function

###################################################################################################3

nms = ["cris.txt", "ever.txt", "grah.txt", "line.txt", "pulc.txt", "sagr.txt", "smar.txt"];
matmatrix = CSV.File(nms, header = 0) |> Tables.matrix
Garray = reshape(matmatrix', 8, 8, 7) ## 8 traits, 7 species
Gvec = [Garray[1:8, 1:8, i] for i in 1:7] ## Vector of Species Gmatrices
DFmats = DataFrame(Gmatrix = Gvec)

xf = DataFrame(XLSX.readtable("Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))

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
display(plot(tree, linewidth=5, tipfont=20))

GancVec = [0.285, 0.118, 0.131, 0.053, 0.212, 0.172, 0.188, 0.210, 0.277,
            0.230, 0.103, 0.074, 0.073, 0.045, 0.023, 0.800, 0.317, 0.054,
            0.121, 0.152, 0.147, 0.654, 0.028, 0.047, 0.124, 0.148, 0.858,
            0.749, 0.581, 0.626, 0.906, 0.553, 0.644, 0.708, 0.731, 0.872]

Ganc =zeros(8, 8)
idx = tril!(trues(size(Ganc)))
Ganc[idx] .= GancVec
Ganc .= Hermitian(Ganc, :L)
uu0 = Ganc #### FullData[7, :Gmatrix]; ## randPosDefMat(8)

alpha = fill(12.0, size(Ganc, 1), size(Ganc, 1));
sigma = fill(1 / sqrt(2) .* 0.10, size(Ganc, 1), size(Ganc, 1));
sigma[diagind(sigma)] .= 0.10;
tspan = (0.0, 1.0)
dt = 0.001
@time test = OUmatrix(uu0, uu0, alpha, sigma, tspan, dt); ## using Ganc as mu and uu0
# Extract the elements from the matrices
test2 = convert(Vector{Matrix{Float64}}, test);
## test2 .= cov2cor.(test2)  ## uncomment to see correlation matrix plots
vals = map(x -> x[tril(trues(size(x)))], test2); 
plot(reuse = false) ## extract lower triangular from test
for i in 1:Integer(dim(test[1], 1) * (dim(test[1], 1) - 1) / 2 + dim(test[1], 1))
plot!(map(x -> x[i], vals), legend = false)
end
current()
