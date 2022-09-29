using DifferentialEquations, Phylo, Plots, Distributions, Statistics,
    Distances, KissABC, JLD2, LinearAlgebra,
    GeometricIntegratorsDiffEq, StochasticIntegrators; pyplot();


        function drift(du, u, p, t)
            alpha, mu, sigma = p
            ## du .= alpha .* u would be BM with drift
            du = alpha * (mu - u) ## OU-like
        end # drift

        function diff(du, u, p, t)
            alpha, mu, sigma = p
            du = sigma ## would be OU
            ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
            ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
      end # diff

    time_tot = 1.0
    tspan = (0.0, time_tot)
u0 = [1.2 0.5;
      0.5 2.1]

alpha1 = 3.0
mu1 = 5.843333 ## Start at the trait means
sigma1 = 1.0
p = (alpha1, mu1, sigma1)


prob = SDEProblem(drift, diff, u0, tspan, p=p)

sol = solve(prob,GIEuler())

npar = 3 ## alpha mu, sigma of the SDE
ndims = length(alpha1)



