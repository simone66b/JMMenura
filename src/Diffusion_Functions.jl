#################################
# Diffusion and Drift functions #
#################################

###################
# Trait functions #
###################

"""
    trait_diff(du, u, p, t)

One possible method for trait drift. Assumes OU process.

Used in SDE problem. (MIGHT NEED TO BE RENAMED)
"""
function trait_drift(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* (mu .- u)
end

"""
    trait_drift_brownian_motion(du, u, p, t)

One possible method for trait drift. Assumes brownian motion.

Used in SDE problem.
"""
function trait_drift_brownian_motion(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* u
end

"""
    trait_diff(du, u, p, t)

One possible method for trait diffusion.

Used in SDE problem.
"""
function trait_diff(du, u, p, t)
    sigma = p.sigma
    du .= sigma
end

"""
    trait_diff_cox_ingersoll_ross_gamma(du, u, p, t)

One possible method for trait diffusion.

Used in SDE problem.
"""
function trait_diff_cox_ingersoll_ross_gamma(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(max.(u,0)) .* sigma
end

"""
    trait_diff_beta(du, u, p, t)

One possible method for trait diffusion.

Used in SDE problem.
"""
function trait_diff_beta(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma
end

####################
# Matrix functions #
####################

"""
    covariance_mat_drift(du, u, p, t)

One possible method for covariance matrix diffusion.

Used in SDE problem.
"""
function matrix_skew_symmetric_drift(du, u, p, t) ## drift function for the SDE
    du = p.a * t * p.A
end

"""
    covariance_mat_diffusion(du, u, p, t)

One possible method for covariance matrix diffusion.

Used in SDE problem.
"""
function matrix_skew_symmetric_diffusion(du, u, p, t) ## diffusion function for the SDE
    for i in 1:size(du)[2]
        du[1:size(du)[1],i] .= p.a .* t.* p.B' 
    end
end


function matrix_OU_drift(du, u, p, t)
    du .= p.alpha .* (p.mu2 - u) ## Mean reversion
end  ## drift function

function matrix_OU_diffusion(du, u, p, t)
    du .= p.sigma ## scaled BM
end