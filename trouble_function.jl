using Phylo, Distributions, Random, Plots, LinearAlgebra, PosDefManifold, DifferentialEquations

function trouble(metric, S, G)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * exp(ℍ(G⁻½ * S * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.expMap function:
              only the Fisher metric is supported for the exponential map"
    end
end

function troubl2(metric, P, G)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * log(ℍ(G⁻½ * P * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.logMap function:
                 only the Fisher metric is supported for the logarithmic map."
    end
end

function matrix_OU_drift(du, u, p, t)
    du .= p.alpha .* (p.mu - u) ## Mean reversion
end  ## drift function

function matrix_OU_diffusion(du, u, p, t)
    du .= p.sigma ## scaled BM
end

n = 9

tspan = (0.0, 0.3)
dt = 0.001

P0 = cor(rand(Wishart(100, Matrix(1I, n, n)  )))

mat_alpha = 12.0 .* ones(n,n)
mat_sigma = (1 / sqrt(2) * 0.10) .* ones(n,n)
mat_mu = copy(P0)

HermId = Hermitian(1.0I(size(P0, 1))) ## Identity matrix
    
## map mean matrix onto tangent space
u0 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(P0), HermId))
mu2 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(mat_mu), HermId))
    
pp = (mu = mu2, alpha = mat_alpha, sigma = mat_sigma) ## tuple of parameters



prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, u0, tspan, p = pp) ## set up the sde problem
sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama

@time temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u);

@time temp = map(x -> exp(x), sol.u);

@time a = map(x -> exp(Hermitian(x)), sol.u);

@time b = exp.(Hermitian.(sol.u));

@time c = [exp(Hermitian(mat)) for mat in sol.u];
