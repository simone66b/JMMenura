using PosDefManifold, LinearAlgebra, DifferentialEquations, StatsBase, Plots
using CSV, Tables, DataFrames, XLSX, Statistics, Pipe
pyplot();

cd("/home/simoneb/Desktop/JMMenura") ## may need to change this to suit
nms = ["cris.txt", "ever.txt", "grah.txt",
"line.txt", "pulc.txt", "sagr.txt", "smar.txt"];
matmatrix = CSV.File(nms, header=0) |> Tables.matrix
Garray = reshape(matmatrix', 8,8,7)
Gvec = [Garray[1:8, 1:8, i] for i in 1:7] ## Species Gmatrices
DFmats = DataFrame(Gmatrix=Gvec)

xf = DataFrame(XLSX.readtable("Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))

vecnames = append!([2], 4:11) # columns to keep
AnolesData = @pipe xf[:,vecnames] |> 
filter(row -> row.Species .== "A. cristatellus" || 
row.Species .== "A. evermanni" ||
row.Species .== "A. grahami" ||
row.Species .== "A. lineatopus" ||
row.Species .== "A. pulchellus" ||
row.Species .== "A. sagrei" ||
row.Species .== "A. smaragdinus", _) |>
groupby(_, :Species) |>
combine(_, 2:9 .=> (x -> mean(skipmissing(x)))) |> ## species means
transform(_, 2:9 .=> (x -> log.(x))) |>
select(_, Not(2:9)) |>
rename(_, names(xf)[vecnames])

FullData = hcat(AnolesData, DFmats)

function OUmatrix(uu0, mu1, alpha, sigma)
     ## starting matrix and mean matrix
function drift(du, u, p, t)
    du .= p.alpha .* (p.mu - u) ## Mean reversion
end  ## drift function

function diffusion(du, u, p, t)
    du .= p.sigma
end

function vec2Hermitian(vec, size)
    Hermitian([ i <= j ? vec[Integer(j * (j - 1)/ 2 + i)] : 0 for i = 1:size, j = 1:size ])
   end

alphavec = alpha[tril!(trues(size(alpha)))]
sigmavec = sigma[tril!(trues(size(sigma)))]

mu1 = Hermitian(mu1) ## mean matrix
uu0 = Hermitian(uu0) ## start matrix
HermId = Hermitian(1.0I(size(uu0, 1))) ## Identity matrix

u0 = logMap(Fisher, uu0, HermId)# map starting matrix onto
                                                       # tangent space (Lie algebra). 
                                                        # must not be Hermitian so we 
                                                        # can mutate the off-diagonals
u1 = u0[tril!(trues(size(u0)))]
mu2 = logMap(Fisher, mu1, HermId) 
mu3 = mu2[tril!(trues(size(mu2)))]
## map mean matrix onto tangent space
pp = (mu = mu3, alpha =  alphavec, sigma = sigmavec) ## tuple of parameters
        
tspan = (0.0, 1.0)
dt = 0.001

prob = SDEProblem(drift, diffusion, u1, tspan, p = pp) # noise_rate_prototype = zeros(2,2)) # set up sde problem
sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama



herms = vec2Hermitian.(sol.u, size(uu0,1))
herms = convert(Vector{Hermitian{Float64, Matrix{Float64}}}, herms)
temp = map(x -> expMap(Fisher, x, HermId), herms) # back transfer trace onto manifold

# replace trace with transformed version.
temp  ## can now plot, print this, etc.
end # OUmatrix function



###################################################################################################3
### FIX THIS 
mu1 = 1000 .* FullData[1,:Gmatrix]

alpha = fill(4.0, size(mu1, 1), size(mu1, 1))
sigma = fill(1/sqrt(2) .* 0.2, size(mu1,1), size(mu1,1))
sigma[diagind(sigma)] .= 0.2
uu0 = randPosDefMat(8)

@time test = OUmatrix(mu1, mu1, alpha, sigma); ## NB using mu1 as starting matrix as well.
# Extract the elements from the matrices

test2 = convert(Vector{Matrix{Float64}}, test)
## test3 = cov2cor.(test2)

plot()
for i in 1:8, j in 1:8
    if i <= j
vals = map(x -> x[i,j], test2)
plot!(vals, legend=false)
    end
end
<<<<<<< HEAD
current()
=======
current()


##################################################3########

############# Data entry ############################3####

using CSV, Tables, DataFrames, XLSX, Statistics, Pipe

cd("/home/simoneb/Desktop/JMMenura") ## may need to change this to suit
nms = ["cris.txt", "ever.txt", "grah.txt",
"line.txt", "pulc.txt", "sagr.txt", "smar.txt"];
matmatrix = CSV.File(nms, header=0) |> Tables.matrix
Garray = reshape(matmatrix', 8,8,7)
Gvec = [Garray[1:8, 1:8, i] for i in 1:7] ## Species Gmatrices
DFmats = DataFrame(Gmatrix=Gvec)

xf = DataFrame(XLSX.readtable("Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli"))

vecnames = append!([2], 4:11) # columns to keep
AnolesData = @pipe xf[:,vecnames] |> 
filter(row -> row.Species .== "A. cristatellus" || 
row.Species .== "A. evermanni" ||
row.Species .== "A. grahami" ||
row.Species .== "A. lineatopus" ||
row.Species .== "A. pulchellus" ||
row.Species .== "A. sagrei" ||
row.Species .== "A. smaragdinus", _) |>
groupby(_, :Species) |>
combine(_, 2:9 .=> (x -> mean(skipmissing(x)))) |> ## species means
transform(_, 2:9 .=> (x -> log.(x))) |>
select(_, Not(2:9)) |>
rename(_, names(xf)[vecnames])

FullData = hcat(AnolesData, DFmats)
>>>>>>> main
