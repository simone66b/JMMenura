using DifferentialEquations, LinearAlgebra, ExponentialUtilities, Statistics;

G = [0 -1; 1 0];
p = (A=G, B=G);
function drift(du, u, p, t) ## drift function for the SDE
    du[1] = t * p.A[1,1]; ## a = 1.0 for example
    du[2] = t * p.A[1,2];
    du[3] = t * p.A[2,1];
    du[4] = t * p.A[2,2];
end # drift

function diffusion(du, u, p, t) ## diffusion function for the SDE
    du[1] = t * p.B[1,1]; ## a = 1.0 for example
    du[2] = t * p.B[1,2];
    du[3] = t * p.B[2,1];
    du[4] = t * p.B[2,2];
end

time_tot = 1.0;
tspan = (0.0, time_tot);
u0 = zeros(2,2);

MM = zeros(4,4)
MM[1,2] = 1
MM[2,1] = -1

fun = SDEFunction(drift, diffusion; mass_matrix=MM);
prob = SDEProblem(fun, diffusion, u0, tspan, p=p); ## setup SDE problem
sol = solve(prob, p=p, dt=0.001); ## Solve using E-M scheme NOT WORKING!
sol.u

zz = last(sol.u)
-zz == zz'

Omega1 = last(sol.u); ## get the final matrix
Omega2 = last(sol.u); ## get another copy

First = exponential!(Omega1); ## matrix exponential
Second = exponential!(-Omega2); ## matrix exponential

result = First * iris * Second; ## reconstruct P_1


##################################################################
using DifferentialEquations, LinearAlgebra, ExponentialUtilities, Statistics;

G = [0 -1; 1 0];
p = (A=G, B=G);

function drift(du, u, p, t) ## drift function for the SDE
    du .= t * p.A;
end # drift

function diffusion(du, u, p, t) ## diffusion function for the SDE
    du .= t * p.B; 
end

time_tot = 1.0;
tspan = (0.0, time_tot);
u0 = zeros(2,2);

prob = SDEProblem(drift, diffusion, u0, tspan, p=p); ## setup SDE problem
sol = solve(prob, ISSEM(theta=1/2, symplectic=true), p=p, dt=0.001);
