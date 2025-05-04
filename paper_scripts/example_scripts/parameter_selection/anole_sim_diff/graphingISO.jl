ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")
priorvec = repeat([Truncated(Normal(0, 50), 0, Inf)], 10)
##priormat = Uniform(0, 2.0)
## push!(priorvec, priormat)
trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix a", "G-matrix b")
###pyplot()
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUISOAnoles50.jld2" pars wts
N = 5000
labels = collect('a':'j')
## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(4, 3, figsize = (15, 12))
subplots_adjust(hspace=0.5)  # Adjust vertical spacing (increase the value for more space)
axes[4, 3].remove()  # Remove the empty subplot
axes[4, 2].remove()  # Remove the empty subplot
for k in 1:10  
    subplot(4,3,k) 
    thisName = nms[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
    ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes, fontsize=12, va="top", ha="right")

parstraitk = [pars[i][k] for i in 1:N]
wtstraitk = [wts[i][k] for i in 1:N]
wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(parstraitk, wtstrait1normalised, 10000, replace=true)

hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins=200)
xlim(0, 10)
title(thisName)
    
if k == 4
    ylabel("Density", fontsize=16)
end

if k > 7
    xlabel(raw"Parameter Value", fontsize=16)
 end
 
 if k == 3
    legend(loc=7)
end
x = range(0.0, 10.0, length=10)
y=pdf(priorvec[k], x)
PyPlot.plot(x, y, color=:red, label="Prior")
end

## show()
PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/ISOanoles50.pdf")