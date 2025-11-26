# Plot Generation Functions for EVA Shiny App
# Contains all plotting logic organized by plot type

library(ggplot2)
library(scales)

# ==================== Plot Theme Configuration ==================== #

# Consistent ggplot2 theme.
get_plot_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
}

# ==================== Metrics Label Helper ==================== #

# Metrics annotation label (n, R², RMSE, MBias, AIC, BIC).
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

# Probability plot vs return period (optional metrics/CI).
create_prob_plot_rperiod <- function(eva_table,
                                     model_rperiods,
                                     model_quant,
                                     distr,
                                     metrics = NULL,
                                     ci_results = NULL,
                                     method = NULL,
                                     ylim = NULL,
                                     data_name = NULL,
                                     title_override = NULL) {
  title_text <- if (!is.null(title_override)) title_override else 
    build_plot_title("Probability Plot", distr, method, data_name)
  p <- ggplot() +
    geom_point(data = eva_table,
               aes(x = rperiod, y = data),
               color = "blue", size = 3, alpha = 0.6) +
    geom_line(aes(x = model_rperiods, y = model_quant),
              color = "red", linewidth = 1.0) +
    scale_x_log10() +
    labs(x = "Return Period",
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

# Probability plot vs exceedance probability.
create_prob_plot_pexc <- function(eva_table,
                                   model_pexc,
                                   model_quant,
                                   distr,
                                   method = NULL,
                                   data_name = NULL,
                                   title_override = NULL) {
  title_text <- if (!is.null(title_override)) title_override else 
    build_plot_title("Probability Plot", distr, method, data_name)
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

# Q-Q plot (sample vs theoretical quantiles).
create_qq_plot <- function(eva_table, distr, params, method = NULL, data_name = NULL, title_override = NULL) {
  title_text <- if (!is.null(title_override)) title_override else 
    build_plot_title("Q-Q Plot", distr, method, data_name)
  
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

# Histogram with fitted PDF overlay.
create_histogram_plot <- function(eva_table, distr, params, method = NULL, data_name = NULL, title_override = NULL) {
  title_text <- if (!is.null(title_override)) title_override else 
    build_plot_title("Histogram with Fitted PDF", distr, method, data_name)
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

# Empirical vs fitted CDF plot.
create_cdf_plot <- function(eva_table, distr, params, method = NULL, data_name = NULL, title_override = NULL) {
  title_text <- if (!is.null(title_override)) title_override else 
    build_plot_title("Empirical vs Fitted CDF", distr, method, data_name)
  
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

# Compare multiple distributions (probability plot).
create_comparison_plot <- function(comparison_results, title_override = NULL) {
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
      x = "Return Period",
      y = "Return Level",
      color = "Distribution"
    ) +
    get_plot_theme() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  if (!is.null(title_override)) {
    p <- p + labs(title = title_override)
  }
  
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

# Compare fitting methods for one distribution.
create_method_comparison_plot <- function(method_comparison_results, title_override = NULL) {
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
      x = "Return Period",
      y = "Return Level",
      color = "Method"
    ) +
    get_plot_theme() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  if (!is.null(title_override)) {
    p <- p + labs(title = title_override)
  }
  
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
