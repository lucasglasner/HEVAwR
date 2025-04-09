source("fit_utils.R")
library("goftest")

# Compute the likelihood value of a given sample and given distribution
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#   log: whether to return the log-likelihood or not. Defaults to FALSE
#
# Returns:
#   A float with the likelihood function value
likelihood <- function(x, distr, params, log = FALSE) {
  dfun <- dprobmodel(x, distr, params) # nolint
  lh <- prod(dfun)
  if (log == TRUE) {
    return(log(lh))
  } else {
    return(lh)
  }
}

# Compute the Akaike Information Criterion (AIC)
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A float with the AIC value
aic_score <- function(x, distr, params) {
  loglh <- likelihood(x, distr, params, log = TRUE)
  return(2 * length(params) - 2 * loglh)
}

# Compute the Bayesian Information Criterion (BIC)
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A float with the AIC value
bic_score <- function(x, distr, params) {
  loglh <- likelihood(x, distr, params, log = TRUE)
  return(log(length(x)) * length(params) - 2 * loglh)
}