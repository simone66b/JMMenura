using DifferentialEquations

##############################
# Matrix and trait evolution #
##############################

"""
    trait_diffusion(x0, tspan, p, mat, trait_drift, trait_diff, dt=0.001) 

Handles evolving a trait for a node
"""
function trait_diffusion(x0, tspan, p, mat, trait_drift, trait_diff, dt=0.001)
    cor1 = cor(mat)
    noise = CorrelatedWienerProcess(cor1, tspan[1],
                                    zeros(dim(cor1)),
                                    zeros(dim(cor1)))
    
    prob = SDEProblem(trait_drift, trait_diff, x0, tspan, p=p, noise=noise);       
    solve(prob, EM(), dt=dt, p=p, adaptive=false)
end

"""
    gen_cov_mat(mat, p, tspan, matrix_drift, u0=zeros(size(mat)), dt = 0.001) 

Handles evolving the covariance matrix of a node
"""
function gen_cov_mat(mat, p, tspan, matrix_drift, u0=zeros(size(mat)), dt = 0.001)
    # Wait why is u0 here. Shouldn't it be mat?
    lowertri = LowerTriangular(mat)
    uppertri = - UpperTriangular(mat)
    skewsymm = lowertri + uppertri
    W = WienerProcess(0.0, 0.0, 0.0)

    pp = (A=skewsymm, B=skewsymm, a=p.a, b=p.b) ## skew symmetric matrices not necessarily the same.
    prob = SDEProblem(matrix_drift, covariance_mat_diffusion, u0, tspan, p=pp, noise=W,
                    noise_rate_prototype=zeros(size(mat))) ## setup SDE problem
    sol = solve(prob, EM(), p=pp, dt=dt)
    Omega1 = exp(last(sol.u)) ## get the final matrix

    Omega1 * mat * Omega1' ## reconstruct P_1
end


function OUmatrix(mat, p, matrix_drift, tspan, dt = 0.001)  
    HermId = Hermitian(1.0I(size(uu0, 1))) ## Identity matrix
    
    ## map mean matrix onto tangent space
    uu0 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(mat), HermId))
    mu2 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(p.mat_mu), HermId))
    
    pp = (mu = mu2, alpha = p.mat_alpha, sigma = p.mat_sigma) ## tuple of parameters
    
    prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, uu0, tspan, p = pp) ## set up the sde problem
    sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
    
    temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
    temp[end] 
    end # OUmatrix function