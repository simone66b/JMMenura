#################################
# Diffusion and Drift functions #
#################################

###################
# Trait functions #
###################

function trait_drift(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* (mu .- u)
end

function trait_drift_brownian_motion(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* u
end

function trait_diff(du, u, p, t)
    sigma = p.sigma
    du .= sigma
end

function trait_diff_cox_ingersoll_ross_gamma(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(u) .* sigma
end

function trait_diff_beta(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma
end

####################
# Matrix functions #
####################

function covariance_mat_drift(du, u, p, t) ## drift function for the SDE
    du .= p.a .* t .* p.A
end

function covariance_mat_diffusion(du, u, p, t) ## diffusion function for the SDE
    du .= p.b .* t .* p.B 
end