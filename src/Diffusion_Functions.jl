#################################
# Diffusion and Drift functions #
#################################

###################
# Trait functions #
###################

"""
    trait_drift_mean_reversion(du, u, p, t)

One possible method for trait drift. Assumes drift of the form

Combined with trait_diffusion_brownian_motion creates the Ornstein-Uhlenbeck process.

Combined with trait_diff_cox_ingersoll_ross_gamma creates the Cox-Ingersoll-Ross process.
"""
function trait_drift_mean_reversion(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* (mu .- u) # The metric should go here Fisher-Rao distance # alpha single scalar 
end

"""
    trait_drift_brownian_motion(du, u, p, t)

One possible method for trait drift. Assumes brownian motion.

Combined with trait_diffusion_brownian_motion creates brownian motion.

Wait is this wrong?
"""
function trait_drift_brownian_motion(du, u, p, t)
    alpha = p.alpha
    mu = p.mu
    du .= alpha .* u
end

"""
    trait_diffusion_brownian_motion(du, u, p, t)

One possible method for trait diffusion.

Combined with trait_drift_mean_reversion creates the Ornstein-Uhlenbeck process.

Combined with trait_drift_brownian_motion creates brownian motion.
"""
function trait_diffusion_brownian_motion(du, u, p, t)
    sigma = p.sigma
    du .= sigma
end

"""
    trait_diff_cox_ingersoll_ross_gamma(du, u, p, t)

One possible method for trait diffusion.

Combined with trait_drift_mean_reversion creates the Cox-Ingersoll-Ross process.
"""
function trait_diffusion_cox_ingersoll_ross_gamma(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(max.(u,0)) .* sigma
end

"""
    trait_diff_beta(du, u, p, t)

One possible method for trait diffusion.

Used in SDE problem.
"""
function trait_diffusion_beta(du, u, p, t)
    sigma = p.sigma;
    du .= sqrt.(abs.(u .* (ones(length(sigma)) .- u))) .* sigma
end

####################
# Matrix functions #
####################

"""
    covariance_mat_drift(du, u, p, t)

One possible method for G matrix drift.

Combined with matrix_diffusion_isospectral evolves the G matrix while keeping the eigenvalues and vectors constant
"""
function matrix_drift_isospectral(du, u, p, t) ## drift function for the SDE
    du = p.a * t * p.A
end

"""
    covariance_mat_diffusion(du, u, p, t)

One possible method for covariance matrix diffusion.

Used in SDE problem.
"""
function matrix_diffusion_isospectral(du, u, p, t) ## diffusion function for the SDE
    for i in 1:size(du)[2]
        du[1:size(du)[1],i] .= p.b .* t.* p.B' 
    end
end


function matrix_drift_mean_reversion(du, u, p, t)
    du .= p.alpha .* (p.mu2 - u) ## Mean reversion
end  ## drift function

function matrix_diffusion_brownian_motion(du, u, p, t)
    du .= p.sigma ## scaled BM
end