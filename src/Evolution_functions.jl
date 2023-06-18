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