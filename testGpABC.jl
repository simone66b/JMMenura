using GpABC
using Distributions
using DifferentialEquations
using Plots
pyplot();
true_params =  [2.0, 1.0, 15.0, 1.0, 1.0, 1.0, 100.0, 1.0, 1.0, 1.0] # nominal parameter values
priors = [Uniform(0.2, 5.), Uniform(0.2, 5.), Uniform(10., 20.),
          Uniform(0.2, 2.), Uniform(0.2, 2.), Uniform(0.2, 2.),
          Uniform(75., 125.), Uniform(0.2, 2.), Uniform(0.2, 2.), 
          Uniform(0.2, 2.)]
param_indices = [1, 3, 9]  #indices of the parameters we want to estimate
priors = priors[param_indices]

#
# ODE solver settings
#
Tspan = (0.0, 10.0)
x0 = [3.0, 2.0, 1.0]
solver = RK4()
saveat = 0.1

GeneReg = function(params::AbstractArray{Float64,1},
    Tspan::Tuple{Float64,Float64}, x0::AbstractArray{Float64,1},
    solver::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm, saveat::Float64)

  if size(params,1) != 10
    throw(ArgumentError("GeneReg needs 10 parameters, $(size(params,1)) were provided"))
  end

  function ODE_3GeneReg(dx, x, par, t)
    dx[1] = par[1]/(1+par[7]*x[3]) - par[4]*x[1]
    dx[2] = par[2]*par[8]*x[1]/(1+par[8]*x[1]) - par[5]*x[2]
    dx[3] = par[3]*par[9]*x[1]*par[10]*x[2]./(1+par[9]*x[1])./(1+par[10]*x[2]) - par[6]*x[3]
  end

  prob = ODEProblem(ODE_3GeneReg, x0 ,Tspan, params)
  Obs = solve(prob, solver, saveat=saveat)

  return hcat(Obs.u...)
end

function simulator_function(var_params)
    params = copy(true_params)
    params[param_indices] .= var_params
    GeneReg(params, Tspan, x0, solver, saveat)
end

reference_data = simulator_function(true_params[param_indices])
plot(reference_data', xlabel="t", ylabel="C(t)", linewidth=2, labels=["u1(t)" "u2(t)" "u3(t)"])

n_particles = 2000
threshold = 1.0
sim_result = SimulatedABCRejection(reference_data, simulator_function, priors, threshold, n_particles;
    max_iter=convert(Int, 2e6),
    write_progress=true)
plot(sim_result)