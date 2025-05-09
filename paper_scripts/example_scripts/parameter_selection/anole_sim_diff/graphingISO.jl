ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")
priorvec = repeat([Truncated(Normal(0, 50), 0, Inf)], 10)
##priormat = Uniform(0, 2.0)
## push!(priorvec, priormat)
trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix a", "G-matrix b")
nms = ["Jaw Length", "Head Width", "Pectoral", "Pelvic", "Humerus", "Ulna", "Femur", "Tibia", "G-matrix a", "G-matrix b"]

###pyplot()
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUISOAnoles50.jld2" pars wts
N = 5000
labels = collect('a':'j')
fig, axes = subplots(5, 2, figsize = (15, 12))
subplots_adjust(hspace=0.5)  # Adjust vertical spacing (increase the value for more space)
##axes[4, 3].remove()  # Remove the empty subplot
##axes[4, 2].remove()  # Remove the empty subplot
for k in 1:10  
    subplot(5,2,k) 
    thisName = nms[k]
    ax = axes[(k-1) ÷ 2 + 1, (k-1) % 2 + 1] 
   ##  ax.text(0.95, 0.95, "($(labels[k]) $thisName)", transform=ax.transAxes, fontsize=12, 
   ##  va="top", ha="right")

parstraitk = [pars[i][k] for i in 1:N]
wtstraitk = [wts[i][k] for i in 1:N]
wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(parstraitk, wtstrait1normalised, 10000, replace=true)

hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins=400)
xlim(0, 10)
title("$(labels[k])) $thisName", loc="left", fontsize=12)
x = range(0.0, 10.0, length=10)
y=pdf(priorvec[k], x)
PyPlot.plot(x, y, color=:red, label="Prior")

    # if iseven(k)
    #     yticks([])
    # end

if k == 5
    ylabel("Density", fontsize=16)
end
if k > 8
    xlabel(raw"Parameter Value", fontsize=16)
else
    xticks([])
 end
 
 if k == 2
    legend(loc=4)
end

end


show()
##PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/ISOanoles50.pdf")