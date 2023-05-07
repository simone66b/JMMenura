#################################
# Diffusion and Drift functions #
#################################


function trait_drift(du, u, p, t)
    alpha = p.alpha;
    mu = p.mu;
    ## du .= alpha .* u would be BM with drift
    du .= alpha .* (mu .- u); ## OU-like
end; # drift

function trait_diff(du, u, p, t)
    sigma = p.sigma;
    du .= sigma; ## would be OU
    ## du .= sqrt.(u) .* sigma ## Cox-Ingersoll-Ross Gamma model
    ## du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma ## Beta model
end; # diff

function covariance_mat_drift(du, u, p, t) ## drift function for the SDE
    du .= p.a .* t .* p.A;
end # drift

function covariance_mat_diffusion(du, u, p, t) ## diffusion function for the SDE
    du .= p.b .* t .* p.B ;
end # diffusion