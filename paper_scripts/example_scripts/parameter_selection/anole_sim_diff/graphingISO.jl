ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")

priorAlpha = Truncated(Normal(0.0, 10.0), 0.0, Inf)
priorvec = repeat([priorAlpha], 10)## 8 traits and one for the matrix_diff
## priorab = Uniform(0, 10)
## priorvec = push!(priorvec, priorab, priorab)

trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix a", "G-Matrix b")
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUISOAnoles.jld2" pars wts
N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i', 'j']

figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(4,3,figsize=(15, 12))

nrow = 4
ncol = 3

for i, j in 1:nrow, 1:ncol
    

for k in 1:10
    thisName = nms[k]
ax = axes[(k-1) ÷ 2 + 1 , (k-1) % 2 + 1] 
parstraitk = [pars[i][k] for i in 1:N]
wtstraitk = [wts[i][k] for i in 1:N]
wtstrait1normalised = Weights(wtstraitk)

samps1 = sample(parstraitk, wtstrait1normalised, 10000, replace=true)
subplot(4,3,k)
hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior")
title(thisName)
ylabel("Density", fontsize=16)
xlabel("Value", fontsize=16)
x = range(0.0, 35.0, length=100)
y=pdf(priorvec[k], x)
PyPlot.plot(x, y, color=:red, label="Prior")

if k == 3
legend(loc=7)
end
end
show()
## PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/ISOAnoles.pdf")