ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")

##priorvec = repeat([Uniform(0, 5)], 8)
## priorvec = push!(priorvec, Uniform(0,10))
##priormat = Uniform(0, 2.0)
## push!(priorvec, priormat)
sigmaPrior = 10.0
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], 9)
trait_data = DataFrame(XLSX.readtable("/home/simoneb/Desktop/JMMenura/anoles_data/Adult measurements for divergence.xlsx", 
"Pmatrix Measurements with outli")) 
nms = push!(names(trait_data)[4:11], "G-Matrix")
###pyplot()
@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUAnoles.jld2" alphas wts
N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
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
hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior",bins=10)
ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes,
            fontsize=12, va="top", ha="right")
            xlim(0, 40)

title(thisName)
    
     if k == 4
        ylabel("Density", fontsize=16)
     end
     if k ==8
        xlabel(raw"α Value", fontsize=16)
     end
    
if k == 9 
    x = range(0.0,30.0, length=100)
else
x = range(0.0, 30.0, length=100)
end
y=pdf(priorvec[k], x)
PyPlot.plot(x, y, color=:red, label="Prior")
 if k == 3
    legend(loc= 7) ## centre right
end
end

PyPlot.savefig("/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/OUanoles.pdf")