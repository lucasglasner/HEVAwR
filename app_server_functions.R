# Server Logic Functions for EVA Shiny App
# Contains core analysis, data processing, and reactive logic

library(shiny)
library(DT)
library(openxlsx)
library(tibble)

# ==================== Data Loading Functions ==================== #

# Load first column from CSV/TXT into numeric vector.
load_data_from_file <- function(filepath) {
  df <- read.csv(filepath, header = FALSE)
  return(df[, 1])
}

# Parse newline-delimited numeric text into vector.
load_data_from_text <- function(text) {
  data <- as.numeric(unlist(strsplit(text, "\n")))
  return(data[!is.na(data)])
}

# ==================== Parameter Control Functions ==================== #

# Return parameter name labels for distribution.
get_param_names <- function(distr) {
  param_names <- switch(distr,
    "norm" = c("Mean", "SD"),
    "lognorm" = c("Mean (log)", "SD (log)"),
    "gamma" = c("Shape", "Rate"),
    "pearson3" = c("Location", "Scale", "Shape"),
    "logpearson3" = c("Location (log)", "Scale (log)", "Shape (log)"),
    "gumbel" = c("Location", "Scale"),
    "gev" = c("Location", "Scale", "Shape")
  )
  return(param_names)
}

# Numeric input controls for fitted parameters.
create_param_controls_ui <- function(fitted_params,
                                      distr) {
  n_params <- length(fitted_params)
  param_names <- get_param_names(distr)
  
  fluidRow(
    lapply(1:n_params, function(i) {
      # First two parameters (location, scale) use 2 decimals
      # Other parameters (shape, etc.) use 3 decimals
      decimals <- if (i <= 2) 2 else 3
      step_size <- if (i <= 2) 0.01 else 0.001
      
      column(4,
        numericInput(paste0("manual_param_", i), 
                    param_names[i], 
                    value = round(fitted_params[i], decimals),
                    step = step_size)
      )
    })
  )
}

# ==================== Analysis Functions ==================== #

# Run EVA analysis; returns model results list.
run_eva_analysis <- function(data,
                              method,
                              distr,
                              target_rperiods,
                              fix_zeros) {
  target_rp <- parse_return_periods(target_rperiods)
  model_rp  <- compute_model_rperiods(target_rp)
  # Run the probability model
  results <- run_probmodel(
    data = data,
    method = method,
    distr = distr,
    model_rperiods = model_rp,
    target_rperiods = target_rp,
    fix_zeros = fix_zeros
  )
  return(results)
}

# Recompute results using manually adjusted params.
recompute_with_manual_params <- function(results,
                                          manual_params,
                                          target_rperiods,
                                          fix_zeros) {
  target_rp <- parse_return_periods(target_rperiods)
  model_rp  <- compute_model_rperiods(target_rp)
  model_pexc <- 1 / model_rp
  # Update parameters
  results$params <- manual_params
  # Recompute model predictions
  model_quant <- get_model_quant(
    model_pexc, results$statistics,
    results$distr, manual_params, fix_zeros
  )
  model_pred <- approx(x = model_rp,
                       y = model_quant,
                       xout = results$eva_table$rperiod)$y
  target_quant <- approx(x = model_rp,
                         y = model_quant,
                         xout = target_rp)$y
  # Recompute metrics
  metrics <- gofmetrics(
    results$eva_table$data, model_pred,
    distr = results$distr, params = manual_params
  )
  # Update results
  model_quant_df <- data.frame(model_quant, row.names = model_rp)
  target_quant_df <- data.frame(target_quant, row.names = target_rp)
  name <- paste(results$distr, results$method, sep = "_")
  colnames(model_quant_df) <- name
  colnames(target_quant_df) <- name
  colnames(metrics) <- name
  results$model_quant <- model_quant_df
  results$target_quant <- target_quant_df
  results$metrics <- metrics
  results$model_pred <- model_pred
  return(results)
}

# ==================== Confidence Interval Functions ==================== #

# Bootstrap confidence intervals for return levels.
compute_bootstrap_ci <- function(data,
                                  n_bootstrap,
                                  distr,
                                  method,
                                  ci_level,
                                  target_rperiods,
                                  statistics,
                                  fix_zeros,
                                  parallel = FALSE) {
  # Run bootstrap
  bootstrap_results <- bootstrap_probmodel(
    x = data,
    niters = n_bootstrap,
    distr = distr,
    method = method,
    parallel = parallel,
    replace = TRUE,
    seed = 123,
    progress = FALSE
  )
  if (is.null(bootstrap_results)) {
    return(NULL)
  }
  # Calculate confidence intervals
  alpha <- 1 - ci_level
  # Get quantiles for model predictions
  target_rp <- parse_return_periods(target_rperiods)
  model_rp  <- compute_model_rperiods(target_rp)
  model_pexc <- 1 / model_rp
  # Compute quantiles for each bootstrap iteration
  n_iters <- nrow(bootstrap_results)
  model_quant_boot <- matrix(
    NA, nrow = length(model_rp), ncol = n_iters
  )
  target_quant_boot <- matrix(
    NA, nrow = length(target_rp), ncol = n_iters
  )
  for (i in 1:n_iters) {
    params_i <- as.numeric(bootstrap_results[i, ])
    if (!any(is.na(params_i))) {
      model_quant_boot[, i] <- get_model_quant(
        model_pexc, statistics, distr, params_i, fix_zeros
      )
      target_quant_boot[, i] <- approx(x = model_rp,
                                       y = model_quant_boot[, i],
                                       xout = target_rp)$y
    }
  }
  # Compute confidence intervals
  model_quant_lower <- apply(
    model_quant_boot, 1, quantile,
    probs = alpha/2, na.rm = TRUE
  )
  model_quant_upper <- apply(
    model_quant_boot, 1, quantile,
    probs = 1 - alpha/2, na.rm = TRUE
  )
  target_quant_lower <- apply(
    target_quant_boot, 1, quantile,
    probs = alpha/2, na.rm = TRUE
  )
  target_quant_upper <- apply(
    target_quant_boot, 1, quantile,
    probs = 1 - alpha/2, na.rm = TRUE
  )
  # Return results
  list(
    model_rp = model_rp,
    model_lower = model_quant_lower,
    model_upper = model_quant_upper,
    target_rp = target_rp,
    target_lower = target_quant_lower,
    target_upper = target_quant_upper,
    ci_level = ci_level
  )
}

# ==================== GOF Test Functions ==================== #

# Run KS, CVM, AD GOF tests; returns summary data frame.
run_gof_tests <- function(data, distr, params) {
  # Run individual GOF tests (Chi-Squared removed)
  ks  <- ks_gof(data, distr, params, alpha = 0.05)
  cvm <- cvm_gof(data, distr, params, alpha = 0.05)
  ad  <- ad_gof(data, distr, params, alpha = 0.05)
  # Build results dataframe
  data.frame(
    Test = c(
      "Kolmogorov-Smirnov",
      "Cramer-von Mises",
      "Anderson-Darling"
    ),
    Statistic = c(
      round(ks$statistic, 4),
      round(cvm$statistic, 4),
      round(ad$statistic, 4)
    ),
    `P-Value` = c(
      round(ks$pvalue, 4),
      round(cvm$pvalue, 4),
      round(ad$pvalue, 4)
    ),
    Passed = c(
      ifelse(ks$test, "✓", "✗"),
      ifelse(cvm$test, "✓", "✗"),
      ifelse(ad$test, "✓", "✗")
    ),
    check.names = FALSE
  )
}

# ==================== Table Functions ==================== #

# Format target quantiles with optional CI columns.
format_quantiles_table <- function(target_quant, ci_results = NULL) {
  quant_df <- target_quant
  # Add confidence intervals if available
  if (!is.null(ci_results)) {
    colnames(quant_df) <- "Estimate"
    quant_df <- cbind(
      quant_df,
      `Lower CI` = ci_results$target_lower,
      `Upper CI` = ci_results$target_upper
    )
  }
  return(quant_df)
}

# ==================== Export Functions ==================== #

# Write Excel report (tables, stats, params, metrics, optional CI/final params).
create_excel_report <- function(results, filename, ci_results = NULL,
                                 initial_params = NULL) {
  # Create workbook
  wb <- createWorkbook()
  # Add worksheets
  addWorksheet(wb, sheetName = "Probability_Table")
  addWorksheet(wb, sheetName = "Statistics")
  addWorksheet(wb, sheetName = "CalibrationParams")
  addWorksheet(wb, sheetName = "FitModel")
  addWorksheet(wb, sheetName = "FitResults")
  addWorksheet(wb, sheetName = "PerformanceMetrics")
  if (!is.null(ci_results)) {
    addWorksheet(wb, sheetName = "ConfidenceIntervals")
  }
  # Prepare parameter data
  param_df <- data.frame(
    distr = paste(results$distr, results$method, sep = "_"),
    loc = ifelse(length(results$params) >= 1, results$params[1], NA),
    scale = ifelse(length(results$params) >= 2, results$params[2], NA),
    shape = ifelse(length(results$params) >= 3, results$params[3], NA)
  )
  # Write data to worksheets
  writeData(wb, sheet = "Probability_Table", x = results$eva_table)
  writeData(wb, sheet = "Statistics", x = results$statistics)
  writeData(wb, sheet = "CalibrationParams", x = param_df)
  # Prepare FitModel data with optional CI
  fitmodel_df <- tibble::rownames_to_column(results$model_quant, "T")
  if (!is.null(ci_results)) {
    fitmodel_df$LowerCI <- ci_results$model_lower
    fitmodel_df$UpperCI <- ci_results$model_upper
  }
  writeData(wb, sheet = "FitModel", x = fitmodel_df)
  writeData(
    wb, sheet = "FitResults",
    x = tibble::rownames_to_column(results$target_quant, "T")
  )
  writeData(
    wb, sheet = "PerformanceMetrics",
    x = tibble::rownames_to_column(results$metrics, "metric")
  )
  # Write confidence intervals if available
  if (!is.null(ci_results)) {
    ci_df <- data.frame(
      ReturnPeriod = ci_results$target_rp,
      Estimate = results$target_quant[, 1],
      LowerCI = ci_results$target_lower,
      UpperCI = ci_results$target_upper,
      CI_Level = ci_results$ci_level
    )
    writeData(wb, sheet = "ConfidenceIntervals", x = ci_df)
  }
  # Write final parameters if they differ from initial
  if (!is.null(initial_params) &&
      !isTRUE(all.equal(initial_params, results$params))) {
    addWorksheet(wb, sheetName = "FinalParams")
    final_param_df <- data.frame(
      distr = paste(results$distr, results$method, sep = "_"),
      loc = ifelse(length(results$params) >= 1, results$params[1], NA),
      scale = ifelse(length(results$params) >= 2, results$params[2], NA),
      shape = ifelse(length(results$params) >= 3, results$params[3], NA)
    )
    writeData(wb, sheet = "FinalParams", x = final_param_df)
  }
  # Save workbook
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}
