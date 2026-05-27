using Pkg
cd("/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated")
ENV["MPLBACKEND"] = "tkagg"  # or "qt5agg"
using PyPlot, JLD2, StatsBase, Distributions, XLSX, DataFrames
PyPlot.matplotlib.use("tkagg")

# Left justify the subtitles to the left with 


@load "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/OU_diff_larger_trait_3/OU_diff_larger_trait_3_importance_result.jld2" tst
alphas = [tst[i][1] for i in 1:5000]
wts = [tst[i][2] for i in 1:5000]

N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
names = ["a) Trait 1", "b) Trait 2", "c) Trait 3", "d) Trait 4", "e) Matrix"]
true_values = [2, 4, 6, 8, 5]

n = 4
sigmaPrior = 50
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], n+1)

## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(2,3,figsize=(15, 8))
for k in 1:5
    thisName = names[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
    alphastraitk = [alphas[i][k] for i in 1:N]
    wtstraitk = [wts[i][k] for i in 1:N]
    wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
    subplot(2,3,k)
    hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins=200)
    # ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes,
    #             fontsize=12, va="top", ha="right")
    title(thisName, loc = "left")
    ax.axvline([true_values[k]], color="orange", label="True value")
    xlim(0, 10)
    
        
    if k == 4
        ylabel("Density", fontsize=16)
    end
    if k ==5
        xlabel("α Value", fontsize=16)
    end
        
    x = range(0, 200, length=300)
    y=pdf(priorvec[k], x)
    PyPlot.plot(x, y, color=:red, label="Prior")
    if k == 3
        legend(loc=1) # PUT IN TOP LEFT
    end
end
axes[6].spines["top"].set_visible(false)
axes[6].spines["bottom"].set_visible(false)
axes[6].spines["left"].set_visible(false) 
axes[6].spines["right"].set_visible(false)
axes[6].xaxis.set_visible(false)
axes[6].yaxis.set_visible(false)

show() ## show the plot
##PyPlot.savefig("OUsim.pdf")


@load "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/BM_model_larger_trait_3/BM_model_larger_trait_3_importance_result.jld2" tst
alphas = [tst[i][1] for i in 1:5000]
wts = [tst[i][2] for i in 1:5000]


N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
names = ["a) Trait 1", "b) Trait 2", "c) Trait 3", "d) Trait 4", "e) Matrix"]
true_values = [2, 4, 6, 8, sqrt(2)]

n = 4
sigmaPrior = 50
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], n+1)

## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(2,3,figsize=(15, 8))
for k in 1:5
    thisName = names[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
    alphastraitk = [alphas[i][k] for i in 1:N]
    wtstraitk = [wts[i][k] for i in 1:N]
    wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
    subplot(2,3,k)
    # ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes,
    #             fontsize=12, va="top", ha="right")
    title(thisName, loc = "left")
    ax.axvline([true_values[k]], color="orange", label="True value")
        
    if k == 4
        ylabel("Density", fontsize=16)
    end
    if k ==5
        xlabel("σ Value", fontsize=16)
    end
        
    if k == 5
        x = range(0,4.0, length=100)
        xlim(0, 4)
        hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior")
    else
        x = range(0, 200, length=300)
        xlim(0, 10)
        hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins=200)
    end
    y=pdf(priorvec[k], x)
    PyPlot.plot(x, y, color=:red, label="Prior")
    if k == 3
        legend(loc=1)
    end
end
axes[6].spines["top"].set_visible(false)
axes[6].spines["bottom"].set_visible(false)
axes[6].spines["left"].set_visible(false) 
axes[6].spines["right"].set_visible(false)
axes[6].xaxis.set_visible(false)
axes[6].yaxis.set_visible(false)

 show()
#PyPlot.savefig("BMsim.pdf")


@load "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/ISO_larger_trait_3/ISO_larger_trait_3_importance_result.jld2" tst
alphas = [tst[i][1] for i in 1:5000]
wts = [tst[i][2] for i in 1:5000]

N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
names = ["a) Trait 1", "b) Trait 2", "c) Trait 3", "d) Trait 4", "e) Matrix a", "f) Matrix b"]
true_values = [2, 4, 6, 8, 5, 5]

n = 4
sigmaPrior = 50
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], n+2)

## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(2,3,figsize=(15, 8))
for k in 1:6
    thisName = names[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
    alphastraitk = [alphas[i][k] for i in 1:N]
    wtstraitk = [wts[i][k] for i in 1:N]
    wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
    subplot(2,3,k)
    hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior", bins = 200)
    # ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes,
    #             fontsize=12, va="top", ha="right")
    title(thisName, loc = "left")
    ax.axvline([true_values[k]], color="orange", label="True value")
    xlim(0, 10)
        
    if k == 4
        ylabel("Density", fontsize=16)
    end
    if k ==5
        xlabel("Value", fontsize=16)
    end
        
    x = range(0, 200, length=300)
    y=pdf(priorvec[k], x)
    PyPlot.plot(x, y, color=:red, label="Prior")
    if k == 3
        legend(loc=1)
    end
end

show()
##PyPlot.savefig("ISOsim.pdf")

@load "/home/simoneb/Desktop/JMMenura/paper_scripts/importance_simulated_scripts/simulated/OU_diff_larger_trait_3/OU_diff_larger_trait_3_importance_result.jld2" tst
alphas = [tst[i][1] for i in 1:5000]
wts = [tst[i][2] for i in 1:5000]

N = 5000
labels = ['a', 'b' ,'c' ,'d' ,'e' ,'f' ,'g' , 'h' , 'i']
names = ["Trait 1", "Trait 2", "Trait 3", "Trait 4", "Matrix"]
true_values = [2, 4, 6, 8, 5]

n = 4
sigmaPrior = 50
prior = Truncated(Normal(0.0, sigmaPrior), 0.0, Inf)
priorvec = repeat([prior], n+1)

## figure(figsize=(30, 20), layout="tight")
fig, axes = subplots(2,3,figsize=(15, 8))
for k in 1:5
    thisName = names[k]
    ax = axes[(k-1) ÷ 3 + 1, (k-1) % 3 + 1] 
    alphastraitk = [alphas[i][k] for i in 1:N]
    wtstraitk = [wts[i][k] for i in 1:N]
    wtstrait1normalised = Weights(wtstraitk)

    samps1 = sample(alphastraitk, wtstrait1normalised, 10000, replace=true)
    subplot(2,3,k)
    hist(samps1; density=true, color="skyblue", edgecolor="black", alpha=0.7, label="Posterior")
    ax.text(0.95, 0.95, "($(labels[k]))", transform=ax.transAxes,
                fontsize=12, va="top", ha="right")
    title(thisName)
    ax.axvline([true_values[k]], color="orange", label="True value")
        
    if k == 4
        ylabel("Density", fontsize=16)
    end
    if k ==5
        xlabel("α Value", fontsize=16)
    end
        
    x = range(0, 200, length=300)
    y=pdf(priorvec[k], x)
    PyPlot.plot(x, y, color=:red, label="Prior")
    if k == 3
        legend(loc=7)
    end
end
axes[6].spines["top"].set_visible(false)
axes[6].spines["bottom"].set_visible(false)
axes[6].spines["left"].set_visible(false) 
axes[6].spines["right"].set_visible(false)
axes[6].xaxis.set_visible(false)
axes[6].yaxis.set_visible(false)

show()
## PyPlot.savefig("OUsim_copy.pdf")