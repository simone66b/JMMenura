ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")
# priorvec = repeat([Uniform(0, 3)], 8) ## for sigma/BM analysis
# priormat = Uniform(0, 2.0) # for sigma/BM analyis

priorvec = repeat([Uniform(0,10)], 9)
## push!(priorvec, priormat)
trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix")
###pyplot()
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/test/AlphaOUAnoles.jld2" alphas wts
N = 5000

## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(3,3,figsize=(15, 12))
for k in 1:9
   
    thisName = nms[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
alphastraitk = [alphas[i][k] for i in 1:N]
wtstraitk = [wts[i][k] for i in 1:N]
wtstrait1normalised = Weights(wtstraitk)

samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
subplot(3,3,k)
hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior")
title(thisName)
    
     if k == 4
        ylabel("Density", fontsize=16)
     end
     if k ==8
        xlabel("Value", fontsize=16)
     end
    
if k == 9 
    x = range(-1.0e-10, 10, length=100)
else
x = range(-1.0e-10, 10, length=100)
end
y=pdf(priorvec[k], x)
plot(x, y, color=:red, label="Prior")
 if k == 3
    legend()
end
end
## show()
PyPlot.savefig("AlphaOUAnoles.pdf")