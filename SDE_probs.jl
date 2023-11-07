using DifferentialEquations, Distributions, Random, LinearAlgebra

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





function drift(du, u, p, t)
    du .= 2*u
end

function diff(du, u, p, t)
    du .= 0.5*u
end

u1 = [1 2; 3 4]
W1 = WienerProcess(0.0, zeros(4))

Random.seed!(1)
prob2 = SDEProblem(drift, diff, u1, tspan, noise = W1, noise_rate_prototype = zeros(2,4))
sol = solve(prob2, EM(), dt = 0.01)


function bruss_f(du, u, p, t)
    du[1] = p[1] + p[2] * u[1] * u[1] * u[2] - p[3] * u[1] - p[4] * u[1]
    du[2] = p[3] * u[1] - p[2] * u[1] * u[1] * u[2]
end
function bruss_g(du, u, p, t)
    du[1, 1] = 0.15 * sqrt(p[1])
    du[1, 2] = 0.15 * sqrt(p[2] * u[1] * u[1] * u[2])
    du[1, 3] = -0.15 * sqrt(p[3] * u[1])
    du[1, 4] = -0.15 * sqrt(p[4] * u[1])
    du[2, 1] = 0
    du[2, 2] = -0.15 * sqrt(p[2] * u[1] * u[1] * u[2])
    du[2, 3] = 0.15 * sqrt(p[3] * 2.5 * u[1])
    du[2, 4] = 0
end
p = (1.0, 1.0, 2.5, 1.0)
"""
Stochastic Brusselator
"""
prob_sde_bruss = SDEProblem(bruss_f, bruss_g, [3.0, 2.0], (0.0, 100.0), p,
    noise_rate_prototype = zeros(2, 4))

solve(prob_sde_bruss)