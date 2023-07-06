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
mu1 = Hermitian(mu1) ## mean matrix
uu0 = Hermitian(uu0) ## start matrix
HermId = Hermitian(1.0I(size(uu0, 1))) ## Identity matrix

u0 = convert(Matrix{Float64}, logMap(Fisher, uu0, HermId))# map starting matrix onto
                                                       # tangent space (Lie algebra). 
                                                        # must not be Hermitian so we 
                                                        # can mutate the off-diagonals

mu = Hermitian(logMap(Fisher, mu1, HermId)) ## map mean matrix onto tangent space
pp = (mu = mu, alpha =  alpha, sigma = sigma) ## tuple of parameters
        
tspan = (0.0, 10.0)
dt = 0.0001

prob = SDEProblem(drift, diffusion, u0, tspan, p=pp) # noise_rate_prototype = zeros(2,2)) # set up sde problem
sol = solve(prob, EM(), p = pp, dt=dt) ## solve using Euler-Maruyama

##temp1 = reshape.(sol.u, 2,2)
temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
sol.u .= temp # replace trace with transformed version.
sol  ## can now plot, print this, etc.
end # OUmatrix function


P0 = [0.329 0.094 -0.083 -0.089 0.293 0.079 0.208 0.268;
0.094 0.449 0.349 0.24 0.071 0.075 0.03 0.009;
-0.083 0.349 1.426 0.487 -0.371 -0.098 -0.053 -0.172;
-0.089 0.24 0.487 0.546 -0.168 0.017 -0.051 -0.081;
0.293 0.071 -0.371 -0.168 1.441 1.008 0.904 0.945;
0.079 0.075 -0.098 0.017 1.008 1.087 0.731 0.78;
0.208 0.03 -0.053 -0.051 0.904 0.731 0.809 0.783;
0.268 0.009 -0.172 -0.081 0.945 0.78 0.783 0.949];

alpha = ones(8,8).*5.0
sigma = ones(8,8)./sqrt(2).*0.25
sigma[diagind(sigma)].=0.25

uu0 = randPosDefMat(8)

@time test = OUmatrix(uu0, P0, alpha, sigma)
# alpha1 = repeat([5.0], 8);
# mu1 = randPosDefMat(8)
# sigma1 = repeat([1.0], 8);



# mu1 = Hermitian([5.0 2.0; 2.0 4.0])
# uu0 =randPosDefMat(2)

# sig = [0.5 0.5/sqrt(2); 0.5/sqrt(2) 0.5]

## test = OUmatrix(uu0, mu1, 5.0, sig);


display(plot(test, legend=false))

test2 = cov2cor.(test.u)
plot(0.0:0001:1.0, test2)
vec1 = map(x -> LowerTriangular(x), test2)
## detvecs = det.(test.u)

display(plot(vec1, reuse=false))

#= function multOffDiag(matrix, multiplier)
    if (ishermitian(matrix))
        matrix = parent(matrix)
    end
    sizeMat = size(matrix, 1)
    result = copy(matrix)
    for i in 1:sizeMat
        for j in 1:sizeMat
            if i != j
                result[i, j] *= multiplier  # Multiply off-diagonal elements by 2
            end
        end     
    end 
    return result
end# multOffDiag  =# 