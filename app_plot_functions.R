# Plot Generation Functions for EVA Shiny App
# Contains all plotting logic organized by plot type

library(ggplot2)
library(scales)

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
                                     ci_results = NULL,
                                     method = NULL,
                                     ylim = NULL,
                                     data_name = NULL) {
  distr_upper <- toupper(distr)
  # Build title with method when provided
  title_base <- if (!is.null(method)) {
    method_label <- switch(method,
      "lmme" = "L-moments",
      "mme" = "Method of Moments",
      "mle" = "Maximum Likelihood",
      method
    )
    paste("Probability Plot -", sprintf("%s (%s)", distr_upper, method_label))
  } else {
    paste("Probability Plot -", distr_upper)
  }
  
  # Prepend data name if provided
  title_text <- if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", title_base)
  } else {
    title_base
  }
  p <- ggplot() +
    geom_point(data = eva_table,
               aes(x = rperiod, y = data),
               color = "blue", size = 3, alpha = 0.6) +
    geom_line(aes(x = model_rperiods, y = model_quant),
              color = "red", linewidth = 1.0) +
    scale_x_log10() +
    labs(x = "Return Period (years)",
         y = "Return Level",
         title = title_text) +
    get_plot_theme()
  
  # Apply custom y-axis limits if provided
  if (!is.null(ylim) && length(ylim) == 2 && !any(is.na(ylim))) {
    p <- p + ylim(ylim[1], ylim[2])
  }
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
#   lang: Language code ("en" or "es"). Default is "en".
# Returns:
#   A ggplot2 object displaying the probability plot with
#   exceedance probability scale and empirical points.
create_prob_plot_pexc <- function(eva_table,
                                   model_pexc,
                                   model_quant,
                                   distr,
                                   method = NULL,
                                   data_name = NULL) {
  distr_upper <- toupper(distr)
  # Build title with method when provided
  title_base <- if (!is.null(method)) {
    method_label <- switch(method,
      "lmme" = "L-moments",
      "mme" = "Method of Moments",
      "mle" = "Maximum Likelihood",
      method
    )
    paste("Probability Plot -", sprintf("%s (%s)", distr_upper, method_label))
  } else {
    paste("Probability Plot -", distr_upper)
  }
  
  # Prepend data name if provided
  title_text <- if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", title_base)
  } else {
    title_base
  }
  ggplot() +
    geom_point(data = eva_table,
               aes(x = pexc, y = data),
               color = "blue", size = 3, alpha = 0.6) +
    geom_line(aes(x = model_pexc, y = model_quant),
             color = "red", linewidth = 1.0) +
    scale_x_log10() +
    labs(x = "Exceedance Probability",
         y = "Return Level",
         title = title_text) +
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
create_qq_plot <- function(eva_table, distr, params, data_name = NULL) {
  distr_upper <- toupper(distr)
  title_base <- paste("Q-Q Plot -", distr_upper)
  title_text <- if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", title_base)
  } else {
    title_base
  }
  
  ggplot(data = eva_table, aes(sample = data)) +
    stat_qq(distribution = function(p) {
      qprobmodel(p, distr, params)
    }, size = 3) +
    stat_qq_line(distribution = function(p) {
      qprobmodel(p, distr, params)
    }, color = "red", linewidth = 1.0) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = title_text) +
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
create_histogram_plot <- function(eva_table, distr, params, data_name = NULL) {
  distr_upper <- toupper(distr)
  title_base <- paste("Histogram with Fitted PDF -", distr_upper)
  title_text <- if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", title_base)
  } else {
    title_base
  }
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
         title = title_text) +
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
create_cdf_plot <- function(eva_table, distr, params, data_name = NULL) {
  distr_upper <- toupper(distr)
  title_base <- paste("Empirical vs Fitted CDF -", distr_upper)
  title_text <- if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", title_base)
  } else {
    title_base
  }
  
  x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
  cdf_vals <- pprobmodel(x_seq, distr, params)
  ggplot() +
    stat_ecdf(data = eva_table, aes(x = data),
             geom = "step", color = "blue", linewidth = 1.0) +
    geom_line(aes(x = x_seq, y = cdf_vals), color = "red", linewidth = 1.0) +
    labs(x = "Value", 
         y = "Cumulative Probability",
         title = title_text) +
    get_plot_theme()
}

# Create a comparison probability plot showing multiple
# distributions.
# Args:
#   comparison_results: Named list where each element is a
#                       distribution result (from
#                       run_eva_analysis).
# Returns:
#   A ggplot object showing empirical data and fitted curves
#   for all distributions.
create_comparison_plot <- function(comparison_results) {
  # Get empirical data from first distribution
  first_result <- comparison_results[[1]]
  eva_table <- first_result$eva_table
  
  # Prepare empirical data
  empirical_data <- data.frame(
    T = eva_table$rperiod,
    Value = eva_table$data
  )
  
  # Prepare model data for all distributions
  model_data_list <- lapply(names(comparison_results), function(distr) {
    result <- comparison_results[[distr]]
    model_table <- result$model_quant
    
    # Extract T values from rownames and Value from column
    T_values <- as.numeric(rownames(model_table))
    Value_values <- model_table[, 1]
    
    data.frame(
      T = T_values,
      Value = Value_values,
      Distribution = toupper(distr)
    )
  })
  
  model_data <- do.call(rbind, model_data_list)
  
  # Define color palette
  n_dists <- length(comparison_results)
  colors <- scales::hue_pal()(n_dists)
  
  # Create base plot
  p <- ggplot() +
    # Model curves
    geom_line(
      data = model_data,
      aes(x = T, y = Value, color = Distribution),
      linewidth = 1.2
    ) +
    # Empirical points
    geom_point(
      data = empirical_data,
      aes(x = T, y = Value),
      size = 2.5,
      alpha = 0.7,
      color = "black"
    ) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(
      x = "Return Period (years)",
      y = "Return Level",
      color = "Distribution"
    ) +
    get_plot_theme() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # ==================== GOF Summary Text Overlay ==================== #
  # Build ✓ / ✗ summary (Chi2, KS, CVM, AD) based on p-value > 0.05.
  has_metrics <- all(vapply(comparison_results, function(r) !is.null(r$metrics), logical(1)))
  if (has_metrics) {
    gof_lines <- c("DIST         KS  CVM AD")
    symbol_for <- function(p) {
      if (is.na(p)) return("NA")
      if (p > 0.05) "✓" else "✗"
    }
    for (distr in names(comparison_results)) {
      mets <- comparison_results[[distr]]$metrics
      pv <- function(name) if (name %in% rownames(mets)) as.numeric(mets[name, 1]) else NA_real_
      line <- sprintf("%-12s %-3s %-3s %-2s", toupper(distr),
              symbol_for(pv("kspvalue")),
                      symbol_for(pv("cvmpvalue")),
                      symbol_for(pv("adpvalue")))
      gof_lines <- c(gof_lines, line)
    }
    overlay_text <- paste(gof_lines, collapse = "\n")
    x_min <- min(model_data$T, na.rm = TRUE)
    y_max <- max(model_data$Value, na.rm = TRUE)
    p <- p + annotate(
      "text", x = x_min, y = y_max, label = overlay_text,
      hjust = 0, vjust = 1, family = "mono", size = 3, lineheight = 0.95,
      fontface = "plain", color = "black"
    )
  }
  
  return(p)
}

# Create a comparison probability plot showing multiple
# fitting methods for a single distribution.
# Args:
#   method_comparison_results: Named list where each element is a
#                              method result (from run_eva_analysis).
# Returns:
#   A ggplot object showing empirical data and fitted curves
#   for all methods.
create_method_comparison_plot <- function(method_comparison_results) {
  # Get empirical data from first method
  first_result <- method_comparison_results[[1]]
  eva_table <- first_result$eva_table
  
  # Prepare empirical data
  empirical_data <- data.frame(
    T = eva_table$rperiod,
    Value = eva_table$data
  )
  
  # Prepare model data for all methods
  model_data_list <- lapply(names(method_comparison_results), function(method) {
    result <- method_comparison_results[[method]]
    model_table <- result$model_quant
    
    # Extract T values from rownames and Value from column
    T_values <- as.numeric(rownames(model_table))
    Value_values <- model_table[, 1]
    
    method_label <- switch(method,
      "lmme" = "L-moments",
      "mme" = "Method of Moments",
      "mle" = "Maximum Likelihood",
      method
    )
    
    data.frame(
      T = T_values,
      Value = Value_values,
      Method = method_label
    )
  })
  
  model_data <- do.call(rbind, model_data_list)
  
  # Create base plot
  p <- ggplot() +
    # Model curves
    geom_line(
      data = model_data,
      aes(x = T, y = Value, color = Method),
      linewidth = 1.2
    ) +
    # Empirical points
    geom_point(
      data = empirical_data,
      aes(x = T, y = Value),
      size = 2.5,
      alpha = 0.7,
      color = "black"
    ) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(
      x = "Return Period (years)",
      y = "Return Level",
      color = "Method"
    ) +
    get_plot_theme() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # ==================== GOF Summary Text Overlay ==================== #
  # Build ✓ / ✗ summary (KS, CVM, AD) based on p-value > 0.05.
  has_metrics <- all(vapply(method_comparison_results, function(r) !is.null(r$metrics), logical(1)))
  if (has_metrics) {
    gof_lines <- c("METHOD  KS  CVM AD")
    symbol_for <- function(p) {
      if (is.na(p)) return("NA")
      if (p > 0.05) "✓" else "✗"
    }
    for (method in names(method_comparison_results)) {
      method_label <- switch(method,
        "lmme" = "L-mom",
        "mme" = "MOM",
        "mle" = "MLE",
        method
      )
      mets <- method_comparison_results[[method]]$metrics
      pv <- function(name) if (name %in% rownames(mets)) as.numeric(mets[name, 1]) else NA_real_
      line <- sprintf("%-7s %-3s %-3s %-2s", method_label,
              symbol_for(pv("kspvalue")),
                      symbol_for(pv("cvmpvalue")),
                      symbol_for(pv("adpvalue")))
      gof_lines <- c(gof_lines, line)
    }
    overlay_text <- paste(gof_lines, collapse = "\n")
    x_min <- min(model_data$T, na.rm = TRUE)
    y_max <- max(model_data$Value, na.rm = TRUE)
    p <- p + annotate(
      "text", x = x_min, y = y_max, label = overlay_text,
      hjust = 0, vjust = 1, family = "mono", size = 3, lineheight = 0.95,
      fontface = "plain", color = "black"
    )
  }
  
  return(p)
}
