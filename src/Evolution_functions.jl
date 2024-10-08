using DifferentialEquations, PosDefManifold

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
                                    zeros(size(cor1)[1]),
                                    zeros(size(cor1)[1]))
    
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


function OUmatrix(mat, para, tspan, matrix_drift, dt = 0.001)
    HermId = Hermitian(1.0I(size(mat, 1))) ## Identity matrix
    
    ## map mean matrix onto tangent space
    uu0 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(mat), HermId))
    mu2 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(para.mat_mu), HermId))
    Ev
    pp = (mu = mu2, alpha = para.mat_alpha, sigma = para.mat_sigma) ## tuple of parameters
    
    # println(matrix_OU_drift)
    # println(matrix_OU_diffusion)
    # println(uu0)
    # println(tspan)
    # println(pp)
    prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, uu0, tspan, p = pp) ## set up the sde problem
    sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
    
    temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
    temp[end] 
    end # OUmatrix function

#########################################
# Functions using each matrix timepoint #
#########################################

function trait_diffusion_each(x0, tspan, p, mats, trait_drift, trait_diff; dt=0.001, small_dt_scale = 10)
    cors1 = cor.(mats)
    u = Vector{Vector{Float64}}()
    t = Vector{Float64}()
    push!(u, x0)
    push!(t, tspan[1])
    for i in 1:(length(cors1)-1)
        cor1 = cors1[i]
        small_tspan = t[end]
        noise = CorrelatedWienerProcess(cor1,small_tspan,
                                    zeros(size(cor1)[1]),
                                    zeros(size(cor1)[1]))
    
        prob = SDEProblem(trait_drift, trait_diff, u[end], (small_tspan, small_tspan + dt), 
                            p=p, noise=noise);       
        sol = solve(prob, EM(), dt=dt/small_dt_scale, p=p, adaptive=false)
        push!(u, sol.u[end])
        push!(t, sol.t[end])
    end
    (u = u, t = t)
end

function OUmatrix_each(mat, para, tspan, matrix_drift, dt = 0.001)
    HermId = Hermitian(1.0I(size(mat, 1))) ## Identity matrix
    
    ## map mean matrix onto tangent space
    uu0 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(mat), HermId))
    mu2 = convert(Matrix{Float64}, logMap(Fisher, Hermitian(para.mat_mu), HermId))
    
    pp = (mu = mu2, alpha = para.mat_alpha, sigma = para.mat_sigma) ## tuple of parameters
    
    prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, uu0, tspan, p = pp) ## set up the sde problem
    sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
    
    temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
    end