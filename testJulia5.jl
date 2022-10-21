using DifferentialEquations, LinearAlgebra, ExponentialUtilities;


## P_0: Classic Fisher's Iris data 
iris =  [0.6856935  -0.0424340    1.2743154   0.5162707
         -0.0424340   0.1899794   -0.3296564  -0.1216394
         1.2743154  -0.3296564    3.1162779   1.2956094
         0.5162707  -0.1216394    1.2956094   0.5810063]

G =  [0.0  0.0424340    -1.2743154   -0.5162707 ## skew symmetric
         -0.0424340   0.0   0.3296564  0.1216394
         1.2743154  -0.3296564    0.0   -1.2956094
         0.5162707  -0.1216394    1.2956094   0.0]


p=(A=G, B=G) ## skew symmetric matrices derived from iris matrix

function drift(du, u, p, t) ## drift function for the SDE
    du .= 1.0 * t .* p.A ## a = 1.0 for example
end # drift

function diffusion(du, u, p, t) ## diffusion function for the SDE
    du .= 1.0 * t .* p.B ## diffusion function b= 1.0 for example
end

time_tot = 1.0;
tspan = (0.0, time_tot);

## u0 = Omega_0
u0 = zeros(4,4);

prob = SDEProblem(drift, diffusion, u0, tspan, p=p); ## setup SDE problem
sol = solve(prob, ISSEM(theta=1/2, symplectic=true), p=p, dt=0.001); ## Solve using symplectic solver

Omega1 = last(sol.u); ## get the final matrix
Omega2neg = -last(sol.u); ## get another copy NB negative

First = exponential!(Omega1); ## matrix exponential
Second = exponential!(Omega2neg); ## matrix exponential

result = First * iris * Second; ## reconstruct P_1
eigen(result).values
eigen(iris).values
result
iris

# eigen(result).values
# 4-element Vector{Float64}:
#  0.023835136854432122
#  0.07820947720632211
#  0.24267075551194184
#  4.228241730427299

# eigen(iris).values
# 4-element Vector{Float64}:
#  0.023835136854432327
#  0.07820947720632193
#  0.24267075551194206
#  4.2282417304273014

# result
# 4×4 Matrix{Float64}:
#   3.00645   -0.917944   0.241816  -1.59811
#  -0.917944   0.451471  -0.14977    0.602206
#   0.241816  -0.14977    0.150938  -0.185023
#  -1.59811    0.602206  -0.185023   0.964095

# iris
# 4×4 Matrix{Float64}:
#   0.685693  -0.042434   1.27432    0.516271
#  -0.042434   0.189979  -0.329656  -0.121639
#   1.27432   -0.329656   3.11628    1.29561
#   0.516271  -0.121639   1.29561    0.581006
