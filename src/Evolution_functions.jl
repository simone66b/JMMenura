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
    ## map mean matrix onto tangent space
    uu0 = convert(Matrix{Float64}, log(Hermitian(mat)))
    mu2 = convert(Matrix{Float64}, log(Hermitian(para.mat_mu)))
    
    pp = (mu2 = mu2, alpha = para.mat_alpha, sigma = para.mat_sigma) ## tuple of parameters

    prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, uu0, tspan, p = pp) ## set up the sde problem
    sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
    
    # temp = map(x -> expMap(Fisher, Hermitian(x), HermId), sol.u) # back transfer trace onto manifold
    # temp[end] 
    Hermitian(exp(Hermitian(sol.u[end])))
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
    ## map mean matrix onto tangent space
    uu0 = convert(Matrix{Float64}, log(Hermitian(mat)))
    mu2 = convert(Matrix{Float64}, log(Hermitian(para.mat_mu)))
    
    pp = (mu2 = mu2, alpha = para.mat_alpha, sigma = para.mat_sigma) ## tuple of parameters
    
    prob = SDEProblem(matrix_OU_drift, matrix_OU_diffusion, uu0, tspan, p = pp) ## set up the sde problem
    sol = solve(prob, EM(), p = pp, dt = dt) ## solve using Euler-Maruyama
    
    temp = exp.(Hermitian.(sol.u)) # back transfer trace onto manifold
end

# function which handles trait evolution 
function trait_evol(;trait_drift = trait_drift_mean_reversion::Function , trait_diffusion = trait_diffusion_brownian_motion::Function, dt = 0.001::Float64, 
                    small_dt_scale = 1.0::Float64)
    function trait_evolving(x0::Vector{Float64}, mat, para::NamedTuple, tspan::Tuple{Float64, Float64}, each::Bool)
        if each 
            cors1 = cor.(mat)
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
            
                prob = SDEProblem(trait_drift, trait_diffusion, u[end], (small_tspan, small_tspan + dt), 
                                    p=para, noise=noise);       
                sol = solve(prob, EM(), dt=dt/small_dt_scale, p=para, adaptive=false)
                push!(u, sol.u[end])
                push!(t, sol.t[end])
            end
            return (u = u, t = t)
        else
            cor1 = (mat[1])
            noise = CorrelatedWienerProcess(cor1, tspan[1],
                                        zeros(size(cor1)[1]),
                                        zeros(size(cor1)[1]))
        
            prob = SDEProblem(trait_drift, trait_diffusion, x0, tspan, p=para, noise=noise);       
            return solve(prob, EM(), dt=dt, p=para, adaptive=false)
        end
    end
end


# This might have to be renamed in future
function mat_evol(;mat_drift = matrix_drift_mean_reversion::Function , mat_diffusion = matrix_diffusion_brownian_motion::Function, dt = 0.001::Float64, mat_err = missing)
    function mat_evolving(mat, para::NamedTuple, tspan::Tuple{Float64, Float64}, each::Bool)
        # println(eigen(mat).values)
        # println()

        # Trying to fix floating point errors
        err = eigen(mat).values[1]
        if err > 0
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat)))
        else
            mat_err_mat = Matrix((min(-10^-12, err))I, size(mat)...)
            # println(eigen(mat - 10*mat_err_mat).values, "\n")
            # println(eigen(para.mu).values)
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat - mat_err_mat)))
        end

        err_mu = eigen(para.mu).values[1]
        if err_mu > 0
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu)))
        else
            mat_err_mat = Matrix((min(-10^-12, err_mu))I, size(para.mu)...)
            # println(eigen(para.mu - 10*mat_err_mat).values, "\n")
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu - mat_err_mat)))
        end
        
        
        para_2 = (mu2 = mu2, para...) # Used to add log mu
    
        prob = SDEProblem(mat_drift, mat_diffusion, uu0, tspan, p = para_2) ## set up the sde problem
        sol = solve(prob, EM(), p = para_2, dt = dt) ## solve using Euler-Maruyama

        timebase = [tspan[1] + 0.005*i for i in 0:ceil((tspan[2] - tspan[1])/ dt)]
        timebase[end] = tspan[2]

        if each
            return (m = exp.(Hermitian.(sol.u)), t = timebase)
        else
            return (m = [exp(Hermitian(sol.u[end]))], t = timebase)
        end        
    end
end


"""
Each must be true for this one. 
"""
function mat_evol_affine(;dt = 0.001::Float64, mat_err = missing, cond_threshold = 1.0e10, verbose = true)
    function mat_evolving(mat, para::NamedTuple, tspan::Tuple{Float64, Float64}, each::Bool)
        err = eigen(mat).values[1]
        if err > 0
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat)))
        else
            mat_err_mat = Matrix((min(-10^-12, err))I, size(mat)...)
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat - 10*mat_err_mat)))
        end

        err_mu = eigen(para.mu).values[1]
        if err_mu > 0
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu)))
        else
            mat_err_mat = Matrix((min(-10^-12, err_mu))I, size(para.mu)...)
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu - 10*mat_err_mat)))
        end
        
        n = size(mat)[1]

        Gs = [Hermitian(Matrix(1.0I, n,n)) for _ in 1:(ceil((tspan[2] - tspan[1])/dt)+1)]

        Gs[1] = Hermitian(mat)
    
    
        for i in 2:length(Gs) 
            W_t = Hermitian(rand(Normal(0,1/sqrt(2)), (n,n)))
            W_t[diagind(W_t)] .*= sqrt(2)
            W_t = Hermitian(W_t)
            last_G = Gs[i-1]
            
            G_cond = cond(last_G)
            if G_cond > cond_threshold
                if verbose
                    @warn "Aborting simulation as condition number $G_cond has exceeded threshold $cond_threshold which 
                    results in instability"
                end
                return (m = nothing, t = nothing), false
            end

            err = real(eigen(last_G).values[1])

            if err > 0
                sqrt_G = Hermitian(sqrt(last_G))
            else
                mat_err_mat = Matrix((-10^-10)I, size(last_G)...)
                sqrt_G = Hermitian(sqrt(last_G - 10*mat_err_mat))
            end

            
            
            inv_sqrt_G = Hermitian(inv(sqrt_G))
            g_mu = inv_sqrt_G*para.mu*inv_sqrt_G
            err = real(eigen(g_mu).values[1])
            if err > 0
                log_g = Hermitian(log(Hermitian(g_mu)))
            else
                mat_err_mat = Matrix((min(-10^-10, err))I, size(g_mu)...)
                log_g = Hermitian(log(Hermitian(g_mu - 10*mat_err_mat)))
            end
            inner = para.alpha*real(log_g)*dt + para.sigma*sqrt(dt)*W_t
            Gs[i] = Hermitian(sqrt_G*exp(Hermitian(inner))*sqrt_G)
        end

        timebase = [tspan[1] + 0.005*i for i in 0:ceil((tspan[2] - tspan[1])/ dt)]
        timebase[end] = tspan[2]

        return (m = Gs, t = timebase), true 
    end
end



function mat_evol_isospectral(;mat_drift = matrix_drift_isospectral::Function , 
                                    mat_diffusion = matrix_diffusion_isospectral::Function, dt = 0.001::Float64, mat_err = missing)
    function mat_evolving_skew_symmetric(mat, para, tspan::Tuple{Float64, Float64}, each::Bool)
        
        n = size(mat)[1]
        u0_proto = zeros(n,n)
        u0 = vcat(u0_proto...)
        
        lowertri = LowerTriangular(mat)
        uppertri = - UpperTriangular(mat)
        skewsymm = lowertri + uppertri
        W = WienerProcess(0.0,zeros(1), 0.0)

        pp = (A=hcat(skewsymm...), B=hcat(skewsymm...), a=para.a, b=para.b) ## skew symmetric matrices not necessarily the same.
        prob = SDEProblem(mat_drift, mat_diffusion, u0, tspan, p=pp, noise=W,
                        noise_rate_prototype=zeros(n^2, 1)) ## setup SDE problem
        sol = solve(prob, EM(), p=pp, dt=dt)

        if each
            res_sol = reshape.(sol.u, n, n)
            Omegas = exp.(res_sol)

            timebase = [tspan[1] + 0.005*i for i in 0:ceil((tspan[2] - tspan[1])/ dt)]
            timebase[end] = tspan[2]

            return (m = [Omega * mat * Omega' for Omega in Omegas], t = timebase)
        else 
            res_sol = reshape(sol.u[end], n, n)
            Omega1 = exp(res_sol) ## get the final matrix

            timebase = [tspan[1] + 0.005*i for i in 0:ceil((tspan[2] - tspan[1])/ dt)]
            timebase[end] = tspan[2]

            return (m = [Omega1 * mat * Omega1'], t = timebase) ## reconstruct P_1
        end
    end
end

"""
Returns two functions. The first is the evolution function for traits which is an empty function while the
second is the function for the matrix.

Each determines if all matrices are stored 
"""
function no_trait_evol(;mat_drift = matrix_drift_mean_reversion::Function , 
    mat_diffusion = matrix_diffusion_brownian_motion::Function, dt = 0.001::Float64, mat_err = missing)

    trait_func = function trait_evolving(x0::Vector{Float64}, mat, para::NamedTuple, tspan::Tuple{Float64, Float64}, each::Bool)
        return (u = nothing, t = nothing)
    end
    mat_func = function mat_evolving(mat, para::NamedTuple, tspan::Tuple{Float64, Float64}, each::Bool)
        # println(eigen(mat).values)
        # println()

        # Trying to fix floating point errors
        err = eigen(mat).values[1]
        if err > 0
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat)))
        else
            mat_err_mat = Matrix((min(-10^-12, err))I, size(mat)...)
            # println(eigen(mat - 10*mat_err_mat).values, "\n")
            # println(eigen(para.mu).values)
            uu0 = convert(Matrix{Float64}, log(Hermitian(mat - mat_err_mat)))
        end

        err_mu = eigen(para.mu).values[1]
        if err_mu > 0
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu)))
        else
            mat_err_mat = Matrix((min(-10^-12, err_mu))I, size(para.mu)...)
            # println(eigen(para.mu - 10*mat_err_mat).values, "\n")
            mu2 = convert(Matrix{Float64}, log(Hermitian(para.mu - mat_err_mat)))
        end
        
        
        para_2 = (mu2 = mu2, para...) # Used to add log mu
    
        prob = SDEProblem(mat_drift, mat_diffusion, uu0, tspan, p = para_2) ## set up the sde problem
        sol = solve(prob, EM(), p = para_2, dt = dt) ## solve using Euler-Maruyama
        if each
            return exp.(Hermitian.(sol.u))
        else
            return [exp(Hermitian(sol.u[end]))]
        end        
    end

    return trait_func, mat_func    
end