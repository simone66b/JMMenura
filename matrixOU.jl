using PosDefManifold, LinearAlgebra, DifferentialEquations, Plots;

function OUmatrix(uu0, mu1, alpha, sigma) ## starting matrix and mean matrix
    function drift(du, u, p, t)
        du .= p.alpha * (p.mu .- u) ## Mean reversion
    end  ## drift function

function diffusion(du, u, p, t)
        du .= p.sigma ## constant sigma
    end # diffusion function

function multOffDiag(matrix, multiplier)
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
end# multOffDiag  

mu1 = Hermitian(mu1) ## mean matrix
uu0 = Hermitian(uu0) ## start matrix
HermId = Hermitian(1.0I(size(uu0)[1])) ## Identity matrix

u0 = convert(Matrix{Float64}, logMap(Fisher, uu0, HermId)) # map starting matrix onto


                                                       # tangent space (Lie algebra). 
                                                        # must not be Hermitian so we 
                                                        # can mutate the off-diagonals

mu = Hermitian(logMap(Fisher, mu1, HermId))
 ## map mean matrix onto tangent space
pp = (mu = mu, alpha =  alpha, sigma = sigma) ## tuple of parameters
        
tspan = (0.0, 10.0)
dt = 0.001
prob = SDEProblem(drift, diffusion, u0, tspan, p=pp) # set up sde problem
sol = solve(prob, EM(), p = pp, dt=dt) ## solve using Euler-Maruyama

temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
sol.u .= temp # replace trace with transformed version.
sol  ## can now plot, print this, etc.
end # OUmatrix function

# uu0=randPosDefMat(5) # random start matrix
# mu1 = randPosDefMat(5)# random mean matrix

mu1 = Hermitian([5.0 2.0; 2.0 4.0])
uu0 = randPosDefMat(2)
@time test = OUmatrix(uu0, mu1, 5.0, 0.1);
plot(test) ## plot the results
