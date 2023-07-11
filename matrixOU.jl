using PosDefManifold, LinearAlgebra, DifferentialEquations, StatsBase, Plots;
pyplot();

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

mu1 = [0.329 0.094 -0.083 -0.089 0.293 0.079 0.208 0.268;
      0.094 0.449 0.349 0.24 0.071 0.075 0.03 0.009;
     -0.083 0.349 1.426 0.487 -0.371 -0.098 -0.053 -0.172;
     -0.089 0.24 0.487 0.546 -0.168 0.017 -0.051 -0.081;
      0.293 0.071 -0.371 -0.168 1.441 1.008 0.904 0.945;
      0.079 0.075 -0.098 0.017 1.008 1.087 0.731 0.78;
      0.208 0.03 -0.053 -0.051 0.904 0.731 0.809 0.783;
      0.268 0.009 -0.172 -0.081 0.945 0.78 0.783 0.949];
      

alpha = fill(0.5, size(mu1, 1), size(mu1, 1))
sigma = fill(1/sqrt(2) .* 0.25, size(mu1,1), size(mu1,1))
sigma[diagind(sigma)] .= 0.25
uu0 = randPosDefMat(8)

@time test = OUmatrix(uu0, mu1, alpha, sigma);
# Extract the elements from the matrices

test2 = convert(Vector{Matrix{Float64}}, test)
test3 = cov2cor.(test2)

plot()
for i in 1:8, j in 1:8
    if i <= j
vals = map(x -> x[i,j], test2)
plot!(vals, legend=false)
    end
end
current()
