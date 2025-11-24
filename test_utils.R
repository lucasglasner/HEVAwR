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


# (Chi-squared GOF test removed from package; retained tests: KS, CVM, AD.)
# Perform a Kolmogorov-Smirnov goodness-of-fit test for the given sample and
# probability model
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#   alpha: Confidence level. Default to 0.05.
#
# Returns:
#   A list with the KS statistic, the pvalue and a bool if the model pass the
#   test or not (0.05 confidence by default).
ks_gof <- function(x, distr, params, alpha = 0.05) {
  ks <- ks.test(x, pprobmodel, distr, params) # nolint
  return(list(
    statistic = ks$statistic,
    pvalue = ks$p.value,
    test = ks$p.value > alpha
  ))
}

# Perform a Cramer-Von-Misses goodness-of-fit test for the given sample and
# probability model
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#   alpha: Confidence level. Default to 0.05.
#
# Returns:
#   A list with the CVM statistic, the pvalue and a bool if the model pass the
#   test or not (0.05 confidence by default).
cvm_gof <- function(x, distr, params, alpha = 0.05) {
  cvm <- cvm.test(x, pprobmodel, distr, params) # nolint
  return(list(
    statistic = cvm$statistic,
    pvalue = cvm$p.value,
    test = cvm$p.value > alpha
  ))
}

# Perform the Anderson-Darling goodness-of-fit test for the given sample and
# probability model
# Args:
#   x: vector with observations
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#   alpha: Confidence level. Default to 0.05.
#
# Returns:
#   A list with the AD statistic, the pvalue and a bool if the model pass the
#   test or not (0.05 confidence by default).
ad_gof <- function(x, distr, params, alpha = 0.05) {
  ad <- ad.test(x, pprobmodel, distr, params) # nolint
  return(list(
    statistic = ad$statistic,
    pvalue = ad$p.value,
    test = ad$p.value > alpha
  ))
}

# Evaluate a probability model in terms of score metrics and goodness of fit
# statistical tests.
# Scores: pearson rsquared, root mean squared error, mean bias, akaike
#         information criterion (AIC), bayesian information criterion (BIC).
# GOF-Tests: Kolmogorov-Smirnov, Cramer-Von Mises and Anderson-Darling.
# Args:
#   x: vector with observations
#   y: vector with simulated values
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#   nbins: Number of intervals to consider when computing the sample histogram.
#          By default it uses the Freedman–Diaconis rule.
#   alpha: Confidence level. Default to 0.05.
#
# Returns:
#   A dataframe with the different scores and statistical test metrics (pvalues
#   and a bool telling if the model is accepted or rejected)
gofmetrics <- function(x, y, distr, params,
                        alpha = 0.05, nbins = NULL) {
  # Deal with missing values
  x <- x[!is.na(x)]
  y <- y[!is.na(x)]
  # Compute score metrics and basic stuff
  n <- length(x)
  r2 <- cor(x, y, method = "pearson") ^ 2
  rmse <- mean((y - x)^2)^0.5
  mbias <- mean(y - x)
  aic <- aic_score(x, distr, params)
  bic <- bic_score(x, distr, params)
  # Goodness-of-fit tests (Chi-squared removed)
  ks <- ks_gof(x, distr, params, alpha = alpha)
  cvm <- cvm_gof(x, distr, params, alpha = alpha)
  ad <- ad_gof(x, distr, params, alpha = alpha)
  # Merge results
  metrics <- c(
    n = n, r2 = r2, rmse = rmse, mbias = mbias,
    aic = aic, bic = bic,
    kspvalue = ks$pvalue,
    cvmpvalue = cvm$pvalue, adpvalue = ad$pvalue,
    kstest = ks$test,
    cvmtest = cvm$test, adtest = ad$test
  )
  metrics <- data.frame(metrics)
  colnames(metrics) <- distr
  return(metrics)
}
