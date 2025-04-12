ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")

priorAlpha = Uniform(0, 3)
priorvec = repeat([priorAlpha], 8)## 8 traits and one for the matrix_diff
priorab = Uniform(0, 10)
priorvec = push!(priorvec, priorab, priorab)

trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix a", "G-Matrix b")
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUISOAnoles.jld2" pars wts
N = 5000

figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(1,2,figsize=(15, 12))
for k in 1:2
    thisName = nms[k+8]
###   global ax = axes[(k-1) ÷ 2 + 1 , (k-1) % 2 + 1] 
global parstraitk = [pars[i][k+8] for i in 1:N]
global wtstraitk = [wts[i][k+8] for i in 1:N]
global wtstrait1normalised = Weights(wtstraitk)

global samps1 = sample(parstraitk, wtstrait1normalised, 10000, replace=true)
subplot(1,2,k)
hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior")
title(thisName)
ylabel("Density", fontsize=16)
xlabel("Value", fontsize=16)
x = range(0.0, 10, length=100)
y=pdf(priorvec[k+8], x)
plot(x, y, color=:red, label="Prior")

if k == 2
legend()
end
end
### show()
PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/ISOAnoles.pdf")