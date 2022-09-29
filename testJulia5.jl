using DifferentialEquations, Phylo, Plots, Distributions, Statistics,
    Distances, KissABC, JLD2, LinearAlgebra, pyplot();


        function drift(du, u, p, t)
            alpha, mu, sigma = p
            ## du .= alpha .* u would be BM with drift
            du .= alpha .* (mu .- u) ## OU-like
        end # drift

        function diff(du, u, p, t)
            alpha, mu, sigma = p
            du .= sigma ## would be OU
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

sol = solve(prob, ImplicitEM(theta=1/2, symplectic=true), p=p, dt=0.001)





   function v(t, q, v, p)
          λ = p[:lambda]
          v[1] = λ*q[1]
          v[2] = λ*q[2]
      end
  
      function B(t, q, B, p, col=0)
          μ = p[:mu]
          if col==0 #whole matrix
              B[1,1] = μ*q[1]
              B[2,1] = μ*q[2]
          elseif col==1
              #just first column
          end
      end
  
      t0 = 0.0
      q0 = [1.0, 1.0]
      lambda  = 2.0
      mu  = 1.0
      p = (lambda=lambda, mu=mu)
  
sde = SDE(1, 1, v, B, t0, q0; parameters=p)


int = Integrator(sde, TableauImplicitEuler(), 0.001)
