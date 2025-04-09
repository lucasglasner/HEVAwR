library("e1071")
library("lmomco")
library("MASS")
library("survival")
library("fitdistrplus")

# ----------------------------- GLOBAL VARIABLES ----------------------------- #
probmodel_names <- c("norm", "lognorm", "gamma", "pearson3",
                     "logpearson3", "gumbel", "gev")
probmodel_nparams <- list(norm = 2, lognorm = 2, gamma = 2, pearson3 = 3,
                          logpearson3 = 3, gumbel = 2, gev = 3)

# ------------------------- Compute empirical moments ------------------------ #
# Compute the first four sample moments of a numeric vector.
# Args: x: A numeric vector containing sample data.
#       type: An integer specifying the type of skewness and kurtosis to
#             compute. Passed to e1071::skewness and e1071::kurtosis. Default
#             is 2 (bias-corrected).
# Returns: A named numeric vector containing mean, variance, skewness, and
#          kurtosis.
moms <- function(x, type = 2) {
  m <- mean(x)
  v <- var(x)
  s <- skewness(x, type = type)
  k <- kurtosis(x, type = type)
  return(c(mean = m, variance = v, skewness = s, kurtosis = k))
}

# Extract a specific sample moment from a numeric vector based on the given
# order.
# Args: x: A numeric vector containing sample data.
#       order: An integer (1 to 4) specifying which moment to return.
#       ...: Additional arguments passed to the moms function.
# Returns: The requested sample moment (mean, variance, skewness, or kurtosis).
emoms <- function(x, order, ...) {
  x_moms <- moms(x, ...)
  if (order == 1) return(x_moms[1])
  if (order == 2) return(x_moms[2])
  if (order == 3) return(x_moms[3])
  if (order == 4) return(x_moms[4])
}

# ---------------------------- Normal distribution --------------------------- #
# Fit a normal distribution to a sample using different methods.
# Args: sample: A numeric vector containing sample data.
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to fitdist (if not using "lmme").
# Returns: A numeric vector containing the estimated parameters
#         (mean, standard deviation).
fit_norm <- function(sample, method = "lmme", ...) {
  if (method == "lmme") {
    params <- parnor(lmoms(sample))$para
  } else {
    params <- fitdist(sample, dist = "norm", method = method, ...)$estimate
  }
  return(as.numeric(params))
}

# -------------------------- Lognormal distribution -------------------------- #
# Fit a log-normal distribution to a sample by applying a normal fit on the
# logarithm of the data.
# Args: sample: A numeric vector containing sample data.
#       ...: Additional arguments passed to fit_norm.
# Returns: A numeric vector containing the estimated parameters (mean,
#          standard deviation) for the log-normal distribution.
fit_lognorm <- function(sample, ...) {
  return(fit_norm(log(sample), ...))
}

# ---------------------------- Gamma distribution ---------------------------- #
# Fit a gamma distribution to a sample using different methods.
# Args: sample: A numeric vector containing sample data.
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to fitdist (if not using "lmme").
# Returns: A numeric vector containing the estimated parameters (shape, rate).
fit_gamma <- function(sample, method = "lmme", ...) {
  if (method == "lmme") {
    params <- pargam(lmoms(sample))$para
    params <- c(params[1], 1 / params[2])
  } else {
    params <- fitdist(sample, dist = "gamma", method = method, ...)$estimate
  }
  return(as.numeric(params))
}

# ----------------------- Pearson type III distribution ---------------------- #
# Compute the probability density function (PDF) of a Pearson Type 3
# distribution.
# Args: x: A numeric vector containing values at which to evaluate the PDF.
#       loc: The location parameter of the Pearson Type 3 distribution.
#       scale: The scale parameter of the Pearson Type 3 distribution.
#       shape: The shape parameter of the Pearson Type 3 distribution.
# Returns: A numeric vector containing the PDF values at the specified points.
dpearson3 <- function(x, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  return(pdfpe3(x, para))
}

# Compute the cumulative distribution function (CDF) of a Pearson Type 3
# distribution.
# Args: q: A numeric vector containing values at which to evaluate the CDF.
#       loc: The location parameter of the Pearson Type 3 distribution.
#       scale: The scale parameter of the Pearson Type 3 distribution.
#       shape: The shape parameter of the Pearson Type 3 distribution.
# Returns: A numeric vector containing the CDF values at the specified points.
ppearson3 <- function(q, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  return(cdfpe3(q, para))
}

# Compute the quantile function (inverse CDF) of a Pearson Type 3 distribution.
# Args: p: A numeric vector of probabilities at which to evaluate the quantiles.
#       loc: The location parameter of the Pearson Type 3 distribution.
#       scale: The scale parameter of the Pearson Type 3 distribution.
#       shape: The shape parameter of the Pearson Type 3 distribution.
# Returns: A numeric vector containing the quantiles for the given
# probabilities.
qpearson3 <- function(p, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  return(quape3(p, para))
}

# Generate random numbers from a Pearson Type 3 distribution.
# Args: n: The number of random values to generate.
#       loc: The location parameter of the Pearson Type 3 distribution.
#       scale: The scale parameter of the Pearson Type 3 distribution.
#       shape: The shape parameter of the Pearson Type 3 distribution.
# Returns: A numeric vector containing `n` random values from the Pearson
# Type 3 distribution.
rpearson3 <- function(n, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  return(rlmomco(n, para))
}

# Compute the statistical moments of a Pearson Type 3 distribution.
# Args: x: A numeric vector of values for which to compute the moments.
#       order: The order of the moment to compute (1 for mean, 2 for variance,
#              etc.).
# Returns: A numeric value corresponding to the computed moment of the
# distribution.
epearson3 <- function(x, order) {
  return(emoms(x, order))
}

# Get the specified moment parameter (mean, scale, or shape) of a Pearson Type
# 3 distribution.
# Args: order: The moment order (1 for mean, 2 for scale, 3 for shape).
#       loc: The location parameter of the Pearson Type 3 distribution.
#       scale: The scale parameter of the Pearson Type 3 distribution.
#       shape: The shape parameter of the Pearson Type 3 distribution.
# Returns: The requested moment parameter.
mpearson3 <- function(order, loc, scale, shape) {
  if (order == 1) {
    return(loc)
  }
  if (order == 2) {
    return(scale)
  }
  if (order == 3) {
    return(shape)
  }
}

# Generate initial parameter estimates for a Pearson Type 3 distribution from
# sample data.
# Args: x: A numeric vector containing sample data.
# Returns: A list containing the estimated location, scale, and shape
# parameters.
fpearson3 <- function(x) {
  loc <- mean(x)
  scale <- sd(x)
  shape <- skewness(x, type = 2)
  fguess <- list(loc = loc, scale = scale, shape = shape)
  return(fguess)
}

# Fit a Pearson Type III distribution to a sample using different methods.
# Args: sample: A numeric vector containing sample data.
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to fitdist (if not using "lmme").
# Returns: A numeric vector containing the estimated parameters (shape, scale,
#          location).
fit_pearson3 <- function(sample, method = "lmme", ...) {
  if (method == "lmme") {
    params <- parpe3(lmoms(sample))$para
  } else if (method == "mme") {
    params <- mmedist(sample, "pearson3", start = fpearson3(sample),
                      order = 1:3, memp = epearson3)$estimate
  } else if (method == "mle") {
    params <- mle2par(sample, type = "pe3")$optim$par
  } else {
    params <- fitdist(sample, dist = "pearson3", method = method,
                      start = fpearson3(sample), ...)$estimate
  }
  return(as.numeric(params))
}

# --------------------- Log pearson type III distribution -------------------- #
# Compute the probability density function (PDF) of a Log-Pearson Type 3
# distribution.
# Args: x: A numeric vector containing values at which to evaluate the PDF.
#       loc: The location parameter of the Log-Pearson Type 3 distribution.
#       scale: The scale parameter of the Log-Pearson Type 3 distribution.
#       shape: The shape parameter of the Log-Pearson Type 3 distribution.
# Returns: A numeric vector containing the PDF values at the specified points.
dlogpearson3 <- function(x, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  pdf  <- pdfpe3(log(x), para) / x
  return(pdf)
}

# Compute the cumulative distribution function (CDF) of a Log-Pearson Type 3
# distribution.
# Args: x: A numeric vector containing values at which to evaluate the CDF.
#       loc: The location parameter of the Log-Pearson Type 3 distribution.
#       scale: The scale parameter of the Log-Pearson Type 3 distribution.
#       shape: The shape parameter of the Log-Pearson Type 3 distribution.
# Returns: A numeric vector containing the CDF values at the specified points.
plogpearson3 <- function(q, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  cdf  <- cdfpe3(log(q), para)
  return(cdf)
}

# Compute the quantile function (inverse CDF) of a Log-Pearson Type 3
# distribution.
# Args: p: A numeric vector of probabilities at which to evaluate the quantiles.
#       loc: The location parameter of the Log-Pearson Type 3 distribution.
#       scale: The scale parameter of the Log-Pearson Type 3 distribution.
#       shape: The shape parameter of the Log-Pearson Type 3 distribution.
# Returns: A numeric vector containing the quantiles for the given
# probabilities.
qlogpearson3 <- function(p, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  quantile <- exp(quape3(p, para))
  return(quantile)
}

# Generate random numbers from a Log-Pearson Type 3 distribution.
# Args: n: The number of random values to generate.
#       loc: The location parameter of the Log-Pearson Type 3 distribution.
#       scale: The scale parameter of the Log-Pearson Type 3 distribution.
#       shape: The shape parameter of the Log-Pearson Type 3 distribution.
# Returns: A numeric vector containing `n` random values from the Log-Pearson
# Type 3 distribution.
rlogpearson3 <- function(n, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "pe3")
  random_x <- exp(rlmomco(n, para))
  return(random_x)
}

# Fit a log-Pearson Type III distribution to a sample by applying a Pearson
# Type III fit on the logarithm of the data.
# Args: sample: A numeric vector containing sample data.
#       ...: Additional arguments passed to fit_pearson3.
# Returns: A numeric vector containing the estimated parameters (shape, scale,
#          location) for the log-Pearson Type III distribution.
fit_logpearson3 <- function(sample, ...) {
  return(fit_pearson3(log(sample), ...))
}


# ---------------------------- Gumbel distribution --------------------------- #
# Compute the probability density function (PDF) of a Gumbel distribution.
# Args: x: A numeric vector containing values at which to evaluate the PDF.
#       loc: The location parameter of the Gumbel distribution.
#       scale: The scale parameter of the Gumbel distribution.
# Returns: A numeric vector containing the PDF values at the specified points.
dgumbel <- function(x, loc, scale) {
  para <- vec2par(c(loc, scale), "gum")
  return(pdfgum(x, para))
}

# Compute the cumulative distribution function (CDF) of a Gumbel distribution.
# Args: q: A numeric vector containing values at which to evaluate the CDF.
#       loc: The location parameter of the Gumbel distribution.
#       scale: The scale parameter of the Gumbel distribution.
# Returns: A numeric vector containing the CDF values at the specified points.
pgumbel <- function(q, loc, scale) {
  para <- vec2par(c(loc, scale), "gum")
  return(cdfgum(q, para))
}

# Compute the quantile function (inverse CDF) of a Gumbel distribution.
# Args: p: A numeric vector of probabilities at which to evaluate the quantiles.
#       loc: The location parameter of the Gumbel distribution.
#       scale: The scale parameter of the Gumbel distribution.
# Returns: A numeric vector containing the quantiles for the given
#          probabilities.
qgumbel <- function(p, loc, scale) {
  para <- vec2par(c(loc, scale), "gum")
  return(quagum(p, para))
}

# Generate random numbers from a Gumbel distribution.
# Args: n: The number of random values to generate.
#       loc: The location parameter of the Gumbel distribution.
#       scale: The scale parameter of the Gumbel distribution.
# Returns: A numeric vector containing `n` random values from the Gumbel
#          distribution.
rgumbel <- function(n, loc, scale) {
  para <- vec2par(c(loc, scale), "gum")
  return(rlmomco(n, para))
}

# Compute the empirical moments of a Gumbel distribution for a given order.
# Args: x: A numeric vector containing sample data.
#       order: The order of the moment to compute (e.g., 1 for the mean, 2
#              for the variance).
# Returns: The value of the specified order moment for the Gumbel distribution
#          based on the empirical data.
egumbel <- function(x, order) {
  return(emoms(x, order))
}

# Compute the theoretical moments of a Gumbel distribution for a given order.
# Args: order: The order of the moment to compute (e.g., 1 for the mean, 2
#              for the variance).
#       loc: The location parameter of the Gumbel distribution.
#       scale: The scale parameter of the Gumbel distribution.
# Returns: The value of the specified order moment for the Gumbel distribution.
mgumbel <- function(order, loc, scale) {
  if (order == 1) {
    m1 <- loc + 0.57721566490153286060 * scale
    return(m1)
  }
  if (order == 2) {
    m2 <- (pi ^ 2 / 6) * (scale ^ 2)
    return(m2)
  }
}

# Generate first guess of Gumbel distribution parameters based on the sample.
# The parameters are estimated using the sample's mean and standard deviation.
# Args: x: A numeric vector containing the sample data.
# Returns: A list with two elements:
#          - loc: The location parameter of the Gumbel distribution.
#          - scale: The scale parameter of the Gumbel distribution.
fgumbel <- function(x) {
  # First guess gumbel parameters from x moments
  scale <- sd(x) / 6 * sqrt(6)
  loc <- mean(x) - 0.57721566490153286060 * scale
  fguess <- list(loc = loc, scale = scale)
  return(fguess)
}

# Fit a Gumbel distribution to a sample using different methods.
# Args: sample: A numeric vector containing sample data.
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to fitdist (if not using "lmme").
# Returns: A numeric vector containing the estimated parameters (location,
#          scale).
fit_gumbel <- function(sample, method = "lmme", ...) {
  if (method == "lmme") {
    params <- pargum(lmoms(sample))$para
  } else if (method == "mme") {
    params <- mmedist(sample, "gumbel", start = fgumbel(sample),
                      order = 1:2, memp = egumbel)$estimate
  } else {
    params <- fitdist(sample, dist = "gumbel", method = method,
                      start = fgumbel(sample), ...)$estimate
  }
  return(as.numeric(params))
}

# ------------------ Generalized extreme value distribution ------------------ #
# Compute the probability density function (PDF) of a Generalized Extreme
# Value (GEV) distribution.
# Args: x: A numeric vector of values for which the PDF should be computed.
#       loc: The location parameter of the GEV distribution.
#       scale: The scale parameter of the GEV distribution.
#       shape: The shape parameter of the GEV distribution.
# Returns: A numeric vector of PDF values corresponding to the provided
#          values based on the specified GEV distribution parameters.
dgev <- function(x, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "gev")
  return(pdfgev(x, para))
}

# Compute the cumulative distribution function (CDF) of a Generalized
# Extreme Value (GEV) distribution.
# Args: q: A numeric vector of quantiles for which the CDF should be computed.
#       loc: The location parameter of the GEV distribution.
#       scale: The scale parameter of the GEV distribution.
#       shape: The shape parameter of the GEV distribution.
# Returns: A numeric vector of CDF values corresponding to the provided
#          quantiles based on the specified GEV distribution parameters.
pgev <- function(q, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "gev")
  return(cdfgev(q, para))
}

# Compute the quantiles of a Generalized Extreme Value (GEV) distribution.
# Args: p: A numeric vector of probabilities for which the quantiles should
#       be computed (values between 0 and 1).
#       loc: The location parameter of the GEV distribution.
#       scale: The scale parameter of the GEV distribution.
#       shape: The shape parameter of the GEV distribution.
# Returns: A numeric vector of quantiles corresponding to the provided
#          probabilities based on the specified GEV distribution parameters.
qgev <- function(p, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "gev")
  return(quagev(p, para))
}

# Generate random samples from a Generalized Extreme Value (GEV) distribution.
# Args: n: The number of random samples to generate.
#       loc: The location parameter of the GEV distribution.
#       scale: The scale parameter of the GEV distribution.
#       shape: The shape parameter of the GEV distribution.
# Returns: A numeric vector of random samples generated from the specified
#          GEV distribution with the given parameters.
rgev <- function(n, loc, scale, shape) {
  para <- vec2par(c(loc, scale, shape), "gev")
  return(rlmomco(n, para))
}

# Compute sample moments for GEV fit routines.
# Args: x: A numeric vector containing the sample data.
#       order: The moment order to compute (1 for mean, 2 for variance,
#              3 for skewness, etc.).
# Returns: The specified moment of the sample, based on the provided order.
egev <- function(x, order) {
  return(emoms(x, order))
}

# Compute the moments of a Generalized Extreme Value (GEV) distribution
# for a given order (1st, 2nd, or 3rd).
# Args: order: The moment order to compute (1 for mean, 2 for variance,
#              3 for skewness).
#       loc: The location parameter of the GEV distribution.
#       scale: The scale parameter of the GEV distribution.
#       shape: The shape parameter of the GEV distribution.
# Returns: The specified moment of the GEV distribution based on the given
#          parameters. Returns Inf if the moment cannot be computed.
mgev <- function(order, loc, scale, shape) {
  shape <- -shape
  if (order == 1) {
    if (shape == 0) {
      return(loc + scale * 0.57721566490153286060)
    } else if (shape < 1) {
      g1 <- gamma(1 - shape)
      return(loc + scale * (g1 - 1) / shape)
    } else {
      return(Inf)
    }
  } else if (order == 2) {
    if (shape == 0) {
      return(scale * (pi ^ 2) / 6)
    } else if (shape < 0.5) {
      g1 <- gamma(1 - shape)
      g2 <- gamma(1 - 2 * shape)
      return((scale ^ 2) * (g2 - g1^2) / (shape ^ 2))
    } else {
      return(Inf)
    }
  } else if (order == 3) {
    if (shape == 0) {
      return(12 * sqrt(6) * 1.202056903159594 / (pi ^ 3))
    } else {
      g1 <- gamma(1 - 1 * shape)
      g2 <- gamma(1 - 2 * shape)
      g3 <- gamma(1 - 3 * shape)
      num <- (g3 - 3 * g2 * g1 + 2 * (g1^3))
      den <- (g2 - (g1 ^ 2)) ^ (3 / 2)
      return(sign(shape) * num / den)
    }
  }
}

# Generate initial parameter guesses for a Generalized Extreme Value (GEV)
# distribution based on L-moments.
# Args: x: A numeric vector containing the sample data.
# Returns: A list containing the initial parameter guesses for the GEV
#          distribution. The list includes:
#          - loc: Location parameter
#          - scale: Scale parameter
#          - shape: Shape parameter
fgev <- function(x) {
  fguess <- pargev(lmoms(x))$par
  fguess <- as.numeric(fguess)
  fguess <- list(loc = fguess[1], scale = fguess[2], shape = fguess[3])
  return(fguess)
}

# Fit a Generalized Extreme Value (GEV) distribution to a sample using
# different methods.
# Args: sample: A numeric vector containing sample data.
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to fitdist (if not using "lmme").
# Returns: A numeric vector containing the estimated parameters (shape,
#          location, scale).
fit_gev <- function(sample, method = "lmme", ...) {
  if (method == "lmme") {
    params <- pargev(lmoms(sample))$para
  } else if (method == "mme") {
    params <- mmedist(sample, "gev", start = fgev(sample),
                      order = 1:3, memp = egev)$estimate
  } else {
    params <- fitdist(sample, dist = "gev", method = method,
                      start = fgev(sample), ...)$estimate
  }
  return(as.numeric(params))
}

# ---------------------------------------------------------------------------- #

# Fit a specified probability distribution to a sample using different methods.
# Args: sample: A numeric vector containing sample data.
#       distr: A character string specifying the distribution. Available options
#              are: "norm", "lognorm", "gamma", "pearson3", "logpearson3",
#              "gumbel", "gev".
#       method: The fitting method to use. Options include:
#               - "lmme" L-moment matching estimation
#               - "mme" Moment matching estimation
#               - "mle" Maximum likelihood estimation
#               - Other methods available in fitdist.
#       ...: Additional arguments passed to the corresponding fitting function.
# Returns: The result of fitting the specified distribution to the sample,
#          using the selected method.
fit_probmodel <- function(sample, distr, method = "lmme", ...) {
  result <- tryCatch({
    # Map distribution names to fitting functions
    distr_map <- list(
      norm        = fit_norm,
      lognorm     = fit_lognorm,
      gamma       = fit_gamma,
      pearson3    = fit_pearson3,
      logpearson3 = fit_logpearson3,
      gumbel      = fit_gumbel,
      gev         = fit_gev
    )
    # Check if the distribution is supported
    if (!distr %in% names(distr_map)) {
      stop(sprintf("Invalid distribution: '%s'. Available options are: %s",
                   distr, paste(names(distr_map), collapse = ", ")))
    }
    # Call the appropriate function
    fit_fun <- distr_map[[distr]]
    return(fit_fun(sample, method = method, ...))
  }, error = function(e) {
    return(rep(NaN, probmodel_nparams[[distr]]))
  })
  return(result)
}

# Calculates the probability density function (PDF) for a given distribution
# and input data.
# Args:
#   x: A numeric vector of values for which the PDF is to be calculated.
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A numeric vector of the calculated PDF values for the input `x` based on the
#   specified distribution and parameters.
dprobmodel <- function(x, distr, params) {
  result <- tryCatch({
    # Map distribution names to density functions
    distr_map <- list(
      norm        = dnorm,
      lognorm     = dlnorm,
      gamma       = dgamma,
      pearson3    = dpearson3,
      logpearson3 = dlogpearson3,
      gumbel      = dgumbel,
      gev         = dgev
    )
    # Check if the distribution is supported
    if (!distr %in% names(distr_map)) {
      stop(sprintf("Invalid distribution: '%s'. Available options are: %s",
                   distr, paste(names(distr_map), collapse = ", ")))
    }
    # Call the appropriate function
    dfun <- distr_map[[distr]]
    args <- c(list(x = x), as.list(params))
    model_density <- do.call(dfun, args)
    return(model_density)
  }, error = function(e) {
    return(rep(NaN, length(x)))
  })
  return(result)
}

# Calculates the cumulative density function (CDF) for a given distribution
# and input data.
# Args:
#   x: A numeric vector of values for which the CDF is to be calculated.
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A numeric vector of the calculated CDF values for the input `x` based on the
#   specified distribution and parameters.
pprobmodel <- function(q, distr, params) {
  result <- tryCatch({
    # Map distribution names to probability functions
    distr_map <- list(
      norm        = pnorm,
      lognorm     = plnorm,
      gamma       = pgamma,
      pearson3    = ppearson3,
      logpearson3 = plogpearson3,
      gumbel      = pgumbel,
      gev         = pgev
    )
    # Check if the distribution is supported
    if (!distr %in% names(distr_map)) {
      stop(sprintf("Invalid distribution: '%s'. Available options are: %s",
                   distr, paste(names(distr_map), collapse = ", ")))
    }
    # Call the appropriate function
    pfun <- distr_map[[distr]]
    args <- c(list(q = q), as.list(params))
    model_probabilities <- do.call(pfun, args)
    return(model_probabilities)
  }, error = function(e) {
    return(rep(NaN, length(q)))
  })
  return(result)
}

# Calculates the inverse cumulative density function or model quantiles for a
# given distribution and input non-excedence probabilities.
# Args:
#   x: A numeric vector of probabilities for quantile calculation
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A numeric vector of the calculated CDF values for the input `x` based on the
#   specified distribution and parameters.
qprobmodel <- function(p, distr, params) {
  result <- tryCatch({
    # Map distribution names to quantile functions
    distr_map <- list(
      norm        = qnorm,
      lognorm     = qlnorm,
      gamma       = qgamma,
      pearson3    = qpearson3,
      logpearson3 = qlogpearson3,
      gumbel      = qgumbel,
      gev         = qgev
    )
    # Check if the distribution is supported
    if (!distr %in% names(distr_map)) {
      stop(sprintf("Invalid distribution: '%s'. Available options are: %s",
                   distr, paste(names(distr_map), collapse = ", ")))
    }
    # Call the appropriate function
    qfun <- distr_map[[distr]]
    args <- c(list(p = p), as.list(params))
    model_quantiles <- do.call(qfun, args)
    return(model_quantiles)
  }, error = function(e) {
    return(rep(NaN, length(p)))
  })
  return(result)
}

# Generate random samples of size n for a given distribution
# Args:
#   n: Random sample size
#   distr: A string specifying the distribution. Must be one of: "norm",
#          "lognorm", "gamma", "pearson3", "logpearson3", "gumbel", or "gev".
#   params: A vector or list of parameters required by the distribution.
#           The number and order of parameters depend on the selected
#           distribution.
#
# Returns:
#   A numeric vector of with the random sample generated with the given
#   distribution
rprobmodel <- function(n, distr, params) {
  result <- tryCatch({
    # Map distribution names to random sample functions
    distr_map <- list(
      norm        = rnorm,
      lognorm     = rlnorm,
      gamma       = rgamma,
      pearson3    = rpearson3,
      logpearson3 = rlogpearson3,
      gumbel      = rgumbel,
      gev         = rgev
    )
    # Check if the distribution is supported
    if (!distr %in% names(distr_map)) {
      stop(sprintf("Invalid distribution: '%s'. Available options are: %s",
                   distr, paste(names(distr_map), collapse = ", ")))
    }
    # Call the appropriate function
    qfun <- distr_map[[distr]]
    args <- c(list(n = n), as.list(params))
    model_random_sample <- do.call(qfun, args)
    return(model_random_sample)
  }, error = function(e) {
    return(rep(NaN, n))
  })
  return(result)
}