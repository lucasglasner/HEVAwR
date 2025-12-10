# Constrained Fitting: Refit parameters with one parameter fixed
# Respects the original fitting method (lmme, mme, or mle)

# Refit all parameters except one (which is held fixed)
# Args:
#   x: Data vector
#   distr: Distribution name
#   method: Fitting method ("lmme", "mme", or "mle")
#   fixed_idx: Index (1-based) of parameter to fix
#   fixed_value: Value to fix the parameter to
# Returns:
#   Vector of refitted parameters with one fixed

refit_with_fixed_param <- function(x, distr, method, fixed_idx, fixed_value) {
  n_params <- probmodel_nparams[[distr]]
  
  if (fixed_idx < 1 || fixed_idx > n_params) {
    stop(sprintf("fixed_idx must be between 1 and %d", n_params))
  }
  
  if (method == "lmme") {
    return(refit_lmom_constrained(x, distr, fixed_idx, fixed_value))
  } else if (method == "mme") {
    return(refit_mom_constrained(x, distr, fixed_idx, fixed_value))
  } else if (method == "mle") {
    return(refit_mle_constrained(x, distr, fixed_idx, fixed_value))
  } else {
    stop("method must be 'lmme', 'mme', or 'mle'")
  }
}

# Refit using L-moments with one parameter fixed
refit_lmom_constrained <- function(x, distr, fixed_idx, fixed_value) {
  n_params <- probmodel_nparams[[distr]]
  free_idx <- setdiff(1:n_params, fixed_idx)
  
  # Special case: if fixing shape to ~0 in GEV, fit as Gumbel for better starting values
  if (distr == "gev" && fixed_idx == 3 && abs(fixed_value) < 0.01) {
    # Fit Gumbel directly to get better starting values
    gumbel_params <- fit_probmodel(x, distr = "gumbel", method = "lmme")
    init_free <- c(gumbel_params[1], gumbel_params[2])
  } else {
    # Get initial parameters from full fit
    init_all <- fit_probmodel(x, distr, method = "lmme")
    init_free <- init_all[free_idx]
  }
  
  # Compute target L-moments from data
  lm_target <- lmoms(x)
  
  # Objective: minimize L-moments distance
  objective <- function(params_free) {
    params_full <- numeric(n_params)
    params_full[free_idx] <- params_free
    params_full[fixed_idx] <- fixed_value
    
    tryCatch({
      # Build para object and compute L-moments
      para_obj <- switch(distr,
        "norm" = list(para = params_full, type = "nor"),
        "lognorm" = list(para = params_full, type = "nor"),
        "gamma" = list(para = c(params_full[1], 1/params_full[2]), type = "gam"),
        "pearson3" = list(para = params_full, type = "pe3"),
        "logpearson3" = list(para = params_full, type = "pe3"),
        "gumbel" = list(para = params_full, type = "gum"),
        "gev" = list(para = params_full, type = "gev"),
        NULL
      )
      
      if (is.null(para_obj)) return(1e10)
      
      class(para_obj) <- para_obj$type
      lm_computed <- par2lmom(para_obj)
      
      # Compare L-moments (use first 3)
      n_lm <- min(3, length(lm_target$lambdas), length(lm_computed$lambdas))
      diff <- sum((lm_target$lambdas[1:n_lm] - lm_computed$lambdas[1:n_lm])^2)
      return(diff)
    }, error = function(e) {
      1e10
    })
  }
  
  # Optimize with multiple random starts to avoid local minima
  best_result <- list(par = init_free, value = Inf)
  
  for (restart in 1:3) {
    # Add some random noise to starting values on restarts
    if (restart == 1) {
      start <- init_free
    } else {
      start <- init_free * (1 + rnorm(length(init_free), 0, 0.1))
    }
    
    opt <- tryCatch({
      optim(start, objective, method = "Nelder-Mead",
            control = list(maxit = 15000, reltol = 1e-12, abstol = 1e-14))
    }, error = function(e) {
      list(par = start, value = Inf, convergence = 1)
    })
    
    if (opt$value < best_result$value) {
      best_result <- opt
    }
  }
  
  # Reconstruct full parameters
  result <- numeric(n_params)
  result[free_idx] <- best_result$par
  result[fixed_idx] <- fixed_value
  return(result)
}

# Refit using ordinary moments with one parameter fixed
refit_mom_constrained <- function(x, distr, fixed_idx, fixed_value) {
  n_params <- probmodel_nparams[[distr]]
  free_idx <- setdiff(1:n_params, fixed_idx)
  
  # Get initial parameters from full fit
  init_all <- fit_probmodel(x, distr, method = "mme")
  init_free <- init_all[free_idx]
  
  # Compute target moments from data
  emp_moms <- moms(x, type = 2)
  
  # Objective: minimize moments distance
  objective <- function(params_free) {
    params_full <- numeric(n_params)
    params_full[free_idx] <- params_free
    params_full[fixed_idx] <- fixed_value
    
    tryCatch({
      # Compute theoretical mean and variance
      theo_mean <- switch(distr,
        "norm" = params_full[1],
        "lognorm" = exp(params_full[1] + params_full[2]^2/2),
        "gamma" = params_full[1] / params_full[2],
        "pearson3" = params_full[1],
        "logpearson3" = params_full[1],
        "gumbel" = params_full[1] + 0.5772*params_full[2],
        "gev" = {
          loc <- params_full[1]
          scale <- params_full[2]
          shape <- params_full[3]
          if (abs(shape) < 1e-6) loc + 0.5772*scale else loc + scale*(gamma(1-shape)-1)/shape
        },
        NA
      )
      
      theo_var <- switch(distr,
        "norm" = params_full[2]^2,
        "lognorm" = (exp(params_full[2]^2) - 1) * exp(2*params_full[1] + params_full[2]^2),
        "gamma" = params_full[1] / params_full[2]^2,
        "pearson3" = params_full[2]^2,
        "logpearson3" = params_full[2]^2,
        "gumbel" = (pi*params_full[2])^2/6,
        "gev" = {
          scale <- params_full[2]
          shape <- params_full[3]
          if (abs(shape) < 1e-6) (pi*scale)^2/6 else scale^2*(gamma(1-2*shape)-gamma(1-shape)^2)/shape^2
        },
        NA
      )
      
      if (is.na(theo_mean) || is.na(theo_var) || theo_var < 0) return(1e10)
      
      diff <- (emp_moms[1] - theo_mean)^2 + (emp_moms[2] - theo_var)^2
      return(diff)
    }, error = function(e) {
      1e10
    })
  }
  
  # Optimize
  opt <- tryCatch({
    optim(init_free, objective, method = "Nelder-Mead",
          control = list(maxit = 10000, reltol = 1e-10, abstol = 1e-12))
  }, error = function(e) {
    list(par = init_free, convergence = 1)
  })
  
  # Reconstruct full parameters
  result <- numeric(n_params)
  result[free_idx] <- opt$par
  result[fixed_idx] <- fixed_value
  return(result)
}

# Refit using MLE with one parameter fixed
refit_mle_constrained <- function(x, distr, fixed_idx, fixed_value) {
  n_params <- probmodel_nparams[[distr]]
  free_idx <- setdiff(1:n_params, fixed_idx)
  
  # Get initial parameters from full fit
  init_all <- fit_probmodel(x, distr, method = "mle")
  init_free <- init_all[free_idx]
  
  # Objective: negative log-likelihood
  objective <- function(params_free) {
    params_full <- numeric(n_params)
    params_full[free_idx] <- params_free
    params_full[fixed_idx] <- fixed_value
    
    # Compute negative log-likelihood
    tryCatch({
      -sum(log(dprobmodel(x, distr, params_full)))
    }, error = function(e) {
      1e10  # Return large value if computation fails
    })
  }
  
  # Optimize using Nelder-Mead
  opt <- tryCatch({
    optim(init_free, objective, method = "Nelder-Mead",
          control = list(maxit = 10000, reltol = 1e-10, abstol = 1e-12))
  }, error = function(e) {
    list(par = init_free, convergence = 1)
  })
  
  # Reconstruct full parameters
  result <- numeric(n_params)
  result[free_idx] <- opt$par
  result[fixed_idx] <- fixed_value
  
  return(result)
}
