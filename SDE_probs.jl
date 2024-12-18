using DifferentialEquations, Distributions, Random, LinearAlgebra, PosDefManifold, GpABC, StatsPlots, ProgressBars, StochasticDiffEq

n = 6

mat = cor(rand(Wishart(100, Matrix(1I, n, n)  )))


u0 = zeros(size(mat))
        
lowertri = LowerTriangular(mat)
uppertri = - UpperTriangular(mat)
skewsymm = lowertri + uppertri
W = WienerProcess(0.0, 0.0, 0.0)

function matrix_skew_symmetric_drift(du, u, p, t) ## drift function for the SDE
    du .= p.a .* t .* p.A
end

function matrix_skew_symmetric_diffusion(du, u, p, t) ## diffusion function for the SDE
    du .= p.b .* t .* p.B 
end
tspan = (0.0, 1.0)

pp = (A=skewsymm, B=skewsymm, a=1, b=1) ## skew symmetric matrices not necessarily the same.
prob = SDEProblem(matrix_skew_symmetric_drift, matrix_skew_symmetric_diffusion, u0, tspan, p=pp, noise=W, noise_rate_prototype=zeros(n,n)) ## setup SDE problem
sol = solve(prob, EM(), p=pp, dt=0.01)





function vector_id(m_sqr, n)
    re_mat = zeros(m_sqr, n)
    for i in 1:n
        re_mat[i+(i-1)*n, i] = 1
    end
    return re_mat
end

function drift(du, u, p, t)
    @show du
    du = p.a * t * p.A
    @show du
    return ones(length(du))
end

function diff(du, u, p, t)
    @show du
    @show t
    @show u
    du = p.a .* t.* p.B 
    @show du
    return du
end

n1 = 3

mat1 = cor(rand(Wishart(100, Matrix(1I, n1, n1)  )))

z1 = zeros(n1, n1)

tspan1 = (0.0, 1.0)
u1 = vcat(mat1...)
W1 = WienerProcess(0.0, ones(1))

lowertri = LowerTriangular(mat1)
uppertri = - UpperTriangular(mat1)
skewsymm = lowertri + uppertri

pp = (A=hcat(skewsymm...), B=hcat(skewsymm...), a=1, b=1)

# Random.seed!(1)
prob2 = SDEProblem(drift, diff, u1, tspan1, pp, noise = W1)
sol = solve(prob2, EM(), dt = 0.01)

res_sol = reshape(sol.u[end], n1, n1)

Omega = exp.(res_sol)

cov_res = Omega * mat1 * Omega'

eig_r = eigvals(cov_res)

eig_s = eigvals(mat1)

eig_r./eig_s



A = zeros(2, 4)
A[1, 1] = 1
A[1, 4] = 1
A[2, 4] = 1

function f(du, u, p, t) 
    @show "drift"
    @show du
    @show size(du)
    @show u
    @show t
    du .= 1.01u
    @show du
end

function g(du, u, p, t)
    @show "diff"
    @show du
    @show size(du)
    @show u
    @show t
    du[1, 1] = 0.3u[1]
    du[1, 4] = 0.12u[2]
    du[2, 4] = 1.8u[2]
    @show du
end

prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = A)
solve(prob, EM(), dt = 0.05)



n1 = 3

mat1 = cor(rand(Wishart(100, Matrix(1I, n1, n1)  )))

lowertri = LowerTriangular(mat1)
uppertri = UpperTriangular(mat1)
skewsymm = lowertri - uppertri
pp = (A=hcat(skewsymm...), B=hcat(skewsymm...), a=1, b=1)

function drift(du, u, p, t)
    # @show "drift"
    # @show du
    # @show u
    # @show p 
    # @show t
    du = p.a * t * p.A
    # @show du
    return ones(length(du))
end

function diff(du, u, p, t)
    # @show "diff"
    # @show du
    # @show size(du)
    # @show u
    # @show p 
    # @show t
    for i in 1:size(du)[2]
        du[1:size(du)[1],i] .= p.a .* t.* p.B' 
    end
    # @show du
    return du
end

z1 = zeros(n1, n1)

tspan1 = (0.0, 1.0)
u1 = vcat(z1...)
W1 = WienerProcess(0.0, zeros(1))

prob2 = SDEProblem(drift, diff, u1, tspan1, pp, noise = W1, noise_rate_prototype = zeros(n1^2, 1))
sol = solve(prob2, EM(), dt = 0.001)

res_sol = reshape(sol.u[end], n1, n1)

Omega = exp(res_sol)

cov_res = Omega * mat1 * Omega'

eig_r = eigvals(cov_res)

eig_s = eigvals(mat1)

eig_r./eig_s