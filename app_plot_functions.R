# Plot Generation Functions for EVA Shiny App
# Contains all plotting logic organized by plot type

library(ggplot2)

# ==================== Plot Theme Configuration ==================== #

# Get a consistent ggplot2 theme for all plots.
# Args:
#   base_size: Numeric value for base font size. Default is 16.
# Returns:
#   A ggplot2 theme object with customized title, axis text, and
#   base font settings.
get_plot_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
}

# ==================== Metrics Label Helper ==================== #

# Create a formatted label displaying goodness-of-fit metrics.
# Args:
#   metrics: Data frame or matrix containing performance metrics
#            (n, R², RMSE, MBias, AIC, BIC).
# Returns:
#   A character string with formatted metric labels for use in
#   plot annotations.
create_metrics_label <- function(metrics) {
  n <- round(metrics["n", 1], 0)
  r2 <- round(metrics["r2", 1], 3)
  rmse <- round(metrics["rmse", 1], 3)
  mbias <- round(metrics["mbias", 1], 3)
  aic <- round(metrics["aic", 1], 2)
  bic <- round(metrics["bic", 1], 2)
  paste0(
    "n = ", n, "\n",
    "R² = ", r2, "\n",
    "RMSE = ", rmse, "\n",
    "MBias = ", mbias, "\n",
    "AIC = ", aic, "\n",
    "BIC = ", bic
  )
}

# ==================== Main Plot Functions ==================== #

# Create a probability plot with return period on the x-axis.
# Args:
#   eva_table: Data frame containing empirical return period and
#              data values.
#   model_rperiods: Numeric vector of model return periods.
#   model_quant: Numeric vector of model-predicted quantiles.
#   distr: Character string specifying the distribution name.
#   metrics: Optional data frame containing goodness-of-fit
#            metrics. Default is NULL.
#   ci_results: Optional list containing confidence interval data.
#               Default is NULL.
# Returns:
#   A ggplot2 object displaying the probability plot with
#   empirical points, fitted curve, optional metrics annotation,
#   and optional confidence bands.
create_prob_plot_rperiod <- function(eva_table,
                                      model_rperiods,
                                      model_quant,
                                      distr,
                                      metrics = NULL,
                                      ci_results = NULL) {
  distr_upper <- toupper(distr)
  p <- ggplot() +
    geom_point(data = eva_table, 
              aes(x = rperiod, y = data), 
              color = "blue", size = 3, alpha = 0.6) +
    geom_line(aes(x = model_rperiods, y = model_quant), 
             color = "red", linewidth = 1.0) +
    scale_x_log10() +
    labs(x = "Return Period (years)", 
         y = "Return Level",
         title = paste("Probability Plot -", distr_upper)) +
    get_plot_theme()
  # Add metrics annotation if provided
  if (!is.null(metrics)) {
    metrics_label <- create_metrics_label(metrics)
    p <- p + annotate(
      "text", x = Inf, y = -Inf,
      label = metrics_label,
      hjust = 1.1, vjust = -0.1,
      size = 5, family = "mono",
      fontface = "plain"
    )
  }
  # Add confidence intervals if available
  if (!is.null(ci_results)) {
    ci_data <- data.frame(
      rperiod = ci_results$model_rp,
      lower = ci_results$model_lower,
      upper = ci_results$model_upper
    )
    p <- p +
      geom_line(
        data = ci_data, aes(x = rperiod, y = lower),
        color = "red", linetype = "dashed",
        linewidth = 0.8
      ) +
      geom_line(
        data = ci_data, aes(x = rperiod, y = upper),
        color = "red", linetype = "dashed",
        linewidth = 0.8
      )
  }
  return(p)
}

# Create a probability plot with exceedance probability on the
# x-axis.
# Args:
#   eva_table: Data frame containing empirical return period and
#              data values.
#   model_pexc: Numeric vector of model exceedance probabilities.
#   model_quant: Numeric vector of model-predicted quantiles.
#   distr: Character string specifying the distribution name.
# Returns:
#   A ggplot2 object displaying the probability plot with
#   exceedance probability scale and empirical points.
create_prob_plot_pexc <- function(eva_table,
                                   model_pexc,
                                   model_quant,
                                   distr) {
  distr_upper <- toupper(distr)
  ggplot() +
    geom_point(data = eva_table, 
              aes(x = pexc, y = data), 
              color = "blue", size = 3, alpha = 0.6) +
    geom_line(aes(x = model_pexc, y = model_quant), 
             color = "red", linewidth = 1.0) +
    scale_x_log10() +
    labs(x = "Exceedance Probability", 
         y = "Return Level",
         title = paste("Probability Plot -", distr_upper)) +
    get_plot_theme()
}

# Create a quantile-quantile (Q-Q) plot comparing sample and
# theoretical quantiles.
# Args:
#   eva_table: Data frame containing sample data values.
#   distr: Character string specifying the distribution name.
#   params: Numeric vector of distribution parameters.
# Returns:
#   A ggplot2 object displaying the Q-Q plot with theoretical
#   quantiles on x-axis and sample quantiles on y-axis.
create_qq_plot <- function(eva_table, distr, params) {
  distr_upper <- toupper(distr)
  ggplot(data = eva_table, aes(sample = data)) +
    stat_qq(distribution = function(p) {
      qprobmodel(p, distr, params)
    }, size = 3) +
    stat_qq_line(distribution = function(p) {
      qprobmodel(p, distr, params)
    }, color = "red", linewidth = 1.0) +
    labs(x = "Theoretical Quantiles", 
         y = "Sample Quantiles",
         title = paste("Q-Q Plot -", distr_upper)) +
    get_plot_theme()
}

# Create a histogram with fitted probability density function
# (PDF) overlay.
# Args:
#   eva_table: Data frame containing sample data values.
#   distr: Character string specifying the distribution name.
#   params: Numeric vector of distribution parameters.
# Returns:
#   A ggplot2 object displaying a histogram with density scale
#   and fitted distribution curve overlay.
create_histogram_plot <- function(eva_table, distr, params) {
  distr_upper <- toupper(distr)
  x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
  pdf_vals <- dprobmodel(x_seq, distr, params)
  # Calculate number of bins using Freedman-Diaconis rule with minimum of 10
  nbins <- max(10, nclass.FD(eva_table$data))
  ggplot() +
    geom_histogram(
      data = eva_table,
      aes(x = data, y = after_stat(density)),
      bins = nbins, fill = "lightblue",
      color = "black", alpha = 0.6,
      linewidth = 0.5
    ) +
    geom_line(aes(x = x_seq, y = pdf_vals), color = "red", linewidth = 1.0) +
    labs(x = "Value", 
         y = "Density",
         title = paste("Histogram with Fitted PDF -", distr_upper)) +
    get_plot_theme()
}

# Create a cumulative distribution function (CDF) plot comparing
# empirical and fitted distributions.
# Args:
#   eva_table: Data frame containing sample data values.
#   distr: Character string specifying the distribution name.
#   params: Numeric vector of distribution parameters.
# Returns:
#   A ggplot2 object displaying the empirical CDF (step function)
#   and fitted theoretical CDF curve.
create_cdf_plot <- function(eva_table, distr, params) {
  distr_upper <- toupper(distr)
  x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
  cdf_vals <- pprobmodel(x_seq, distr, params)
  ggplot() +
    stat_ecdf(data = eva_table, aes(x = data), 
             geom = "step", color = "blue", linewidth = 1.0) +
    geom_line(aes(x = x_seq, y = cdf_vals), color = "red", linewidth = 1.0) +
    labs(x = "Value", 
         y = "Cumulative Probability",
         title = paste("Empirical vs Fitted CDF -", distr_upper)) +
    get_plot_theme()
}
