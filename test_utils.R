source("fit_utils.R")
library("goftest")

# Likelihood (optionally log) for sample under distribution.
likelihood <- function(x, distr, params, log = FALSE) {
  dfun <- dprobmodel(x, distr, params) # nolint
  lh <- prod(dfun)
  if (log == TRUE) {
    return(log(lh))
  } else {
    return(lh)
  }
}

# AIC score = 2k - 2 logL.
aic_score <- function(x, distr, params) {
  loglh <- likelihood(x, distr, params, log = TRUE)
  return(2 * length(params) - 2 * loglh)
}

# BIC score = k log(n) - 2 logL.
bic_score <- function(x, distr, params) {
  loglh <- likelihood(x, distr, params, log = TRUE)
  return(log(length(x)) * length(params) - 2 * loglh)
}


# (Chi-squared GOF test removed from package; retained tests: KS, CVM, AD.)
# Kolmogorov-Smirnov GOF test wrapper.
ks_gof <- function(x, distr, params, alpha = 0.05) {
  ks <- suppressWarnings(ks.test(x, pprobmodel, distr, params)) # nolint
  return(list(
    statistic = ks$statistic,
    pvalue = ks$p.value,
    test = ks$p.value > alpha
  ))
}

# Cramer–Von Mises GOF test wrapper.
cvm_gof <- function(x, distr, params, alpha = 0.05) {
  cvm <- cvm.test(x, pprobmodel, distr, params) # nolint
  return(list(
    statistic = cvm$statistic,
    pvalue = cvm$p.value,
    test = cvm$p.value > alpha
  ))
}

# Anderson-Darling GOF test wrapper.
ad_gof <- function(x, distr, params, alpha = 0.05) {
  ad <- ad.test(x, pprobmodel, distr, params) # nolint
  return(list(
    statistic = ad$statistic,
    pvalue = ad$p.value,
    test = ad$p.value > alpha
  ))
}

# Compute metrics (r2, rmse, bias, AIC, BIC) + GOF (KS, CVM, AD).
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
