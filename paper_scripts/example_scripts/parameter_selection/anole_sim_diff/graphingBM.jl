ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")
priorvec = repeat([Truncated(Normal(0, 50), 0, Inf)], 9)
##priormat = Uniform(0, 2.0)
## push!(priorvec, priormat)
trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix")
###pyplot()
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/BManoles50.jld2" sigmas wts
N = 5000
## labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
labels = collect('a':'i')
## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(3,3,figsize=(15, 12))
for k in 1:9
   
    thisName = nms[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
alphastraitk = [sigmas[i][k] for i in 1:N]
wtstraitk = [wts[i][k] for i in 1:N]
wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
    subplot(3,3,k)
    if k == 9
        bins=10
        xlim(0, 2)
        x = range(0.0, 2.0, length=10)
else
    bins=200
    xlim(0, 10)
    x = range(0.0, 10.0, length=10)
end
hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins=bins)
ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes, fontsize=12, va="top", ha="right")

title(thisName)
    
if k == 4
    ylabel("Density", fontsize=16)
end

if k ==8
    xlabel(raw"σ Value", fontsize=16)
 end
    
y=pdf(priorvec[k], x)
PyPlot.plot(x, y, color=:red, label="Prior")
 if k == 3
    legend(loc=7)
end
end

### show()
PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/BManoles50.pdf")