source("fit_utils.R")

# ==================== Global Helper Utilities (Refactored) ==================== #

# Human-readable method label
get_method_label <- function(method) {
  switch(method,
    "lmme" = "L-moments",
    "mme"  = "Method of Moments",
    "mle"  = "Maximum Likelihood",
    method
  )
}

# Sanitize data name for filenames
sanitize_data_name <- function(name, default = "EVA_Data") {
  if (is.null(name) || identical(name, "")) return(default)
  gsub("[^A-Za-z0-9_-]", "_", name)
}

# Parse comma-separated return periods into numeric vector
parse_return_periods <- function(rp_str) {
  as.numeric(unlist(strsplit(rp_str, ",")))
}

# Generate expanded model return periods for smooth plotting
compute_model_rperiods <- function(target_rp, length_out = 100) {
  max_target <- max(target_rp, na.rm = TRUE)
  exp(seq(log(1.01), log(max_target * 2), length.out = length_out))
}

# Build standardized plot title: Data | Kind - DIST (Method)
build_plot_title <- function(plot_kind, distr, method = NULL, data_name = NULL) {
  distr_upper <- toupper(distr)
  base <- if (!is.null(method)) {
    paste0(plot_kind, " - ", distr_upper, " (", get_method_label(method), ")")
  } else {
    paste0(plot_kind, " - ", distr_upper)
  }
  if (!is.null(data_name) && data_name != "") {
    paste(data_name, "|", base)
  } else base
}



# Extract non-NA values and their corresponding IDs from a
# vector or data frame column.
# Args: data: A numeric vector or data frame containing sample data.
#       col: (Optional) The column name to extract from if data is a data frame.
# Returns: A list with two elements:
#          - ndata: The non-NA values.
#          - ids: The corresponding row names or indices of the non-NA values.
extract_non_na <- function(data, col = NULL) {
  if (is.data.frame(data)) {
    # If col is not provided and data has only one column, use that column
    if (is.null(col)) {
      if (ncol(data) != 1) {
        stop(
          "Column name must be provided for data",
          "with multiple columns."
        )
      }
      col <- colnames(data)[1]
    }
    ndata <- data[[col]]
    ndata <- ndata[!is.na(ndata)]
    ids <- rownames(data)[!is.na(data)]
  } else {
    # If data is a vector
    ndata <- data[!is.na(data)]
    ids <- seq_along(data)[!is.na(data)]
  }
  list(ndata = ndata, ids = ids)
}


# Compute sample statistics (mean, sd, skewness, kurtosis, min, max,
# zero probability) for a vector or data frame column.
# Args: data: A numeric vector or data frame containing sample data.
#       col: (Optional) The column name to extract from if data is a data frame.
#       remove_zeros: Logical; if TRUE, zero values are removed before
#                     computing statistics.
# Returns: A matrix with named statistics as columns.
get_sample_statistics <- function(data,
                                   col = NULL,
                                   remove_zeros = FALSE) {
  extracted <- extract_non_na(data, col)
  ndata <- extracted$ndata
  # Calculate statistics
  pnonzero <- 1 - sum(ndata == 0) / length(ndata)
  if (remove_zeros) {
    ndata <- ndata[ndata > 0]
  }
  sample_moms <- as.list(moms(ndata))
  ndata_mean <- sample_moms$mean
  ndata_sd <- sqrt(sample_moms$variance)
  ndata_skw <- sample_moms$skewness
  ndata_kur <- sample_moms$kurtosis
  ndata_min <- min(ndata)
  ndata_max <- max(ndata)
  stats <- list(mean = ndata_mean, sd = ndata_sd,
                skew = ndata_skw, kurt = ndata_kur,
                min = ndata_min, max = ndata_max,
                pnonzero = pnonzero)
  return(as.data.frame(stats))
}

# Build an exceedance probability table for a vector or data frame column.
# Args: data: A numeric vector or data frame containing sample data.
#       col: (Optional) The column name to extract from if data is a data frame.
#       remove_zeros: Logical; if TRUE, zero values are removed before building
#                     the table.
# Returns: A data frame with columns: rank, data, id, pexc (empirical
#          probability), rperiod (recurrence period).
build_eva_table <- function(data,
                             col = NULL,
                             remove_zeros = FALSE) {
  extracted <- extract_non_na(data, col)
  ndata <- extracted$ndata
  ids <- extracted$ids
  # Remove zero values
  pnonzero <- 1 - sum(ndata == 0) / length(ndata)
  if (remove_zeros) {
    ndata <- ndata[ndata > 0]
    ids <- ids[ndata > 0]
  }
  norder <- order(ndata, decreasing = TRUE) # Order the data
  ids  <- ids[norder]
  ndata  <- ndata[norder]
  # Calculate ranks and exceedance probabilities
  ndata_rank <- rank(-ndata, ties.method = "max")
  ndata_pexc <- ndata_rank / (length(ndata) + 1) # Empirical probability
  if (pnonzero < 1) {
    ndata_pexc <- ndata_pexc * pnonzero # Adjust for zero values
  }

  eva_table <- list(rank = ndata_rank, data = ndata,
                    id = ids, pexc = ndata_pexc, rperiod = 1 / ndata_pexc)
  eva_table <- as.data.frame(eva_table)
  return(eva_table)
}


# Compute model quantiles for given exceedance probabilities.
# Args:
#   model_pexc: Numeric vector of exceedance probabilities for the model.
#   model_rperiods: Numeric vector of recurrence periods for the model.
#   statistics: Data frame or list of sample statistics (must include pnonzero).
#   distr: Character string specifying the distribution name.
#   params: List or vector of fitted distribution parameters.
#   fix_zeros: Logical; if TRUE, adjusts quantiles for zero-inflated data.
# Returns:
#   Numeric vector of model quantiles corresponding to model_rperiods.
get_model_quant <- function(model_pexc,
                             statistics,
                             distr,
                             params,
                             fix_zeros) {
  model_rperiods <- 1 / model_pexc
  model_quant <- qprobmodel(
    (1 - model_pexc), distr = distr, params = params
  )
  if (fix_zeros) {
    model_quant <- approx(x = model_rperiods / statistics$pnonzero,
                          y = model_quant,
                          xout = model_rperiods)$y
  }
  model_quant[model_quant <= 0]   <- 0 # Ensure no negative quantiles
  model_quant[is.na(model_quant)] <- 0 # Handle NA values
  return(model_quant)
}


# Fit a probability model to exceedance data and compute quantiles and metrics.
# Args:
#   data: Numeric vector or data frame containing sample data.
#   method: Character string specifying the fitting method.
#   distr: Character string specifying the distribution name.
#   model_rperiods: Numeric vector of recurrence periods for the model.
#   target_rperiods: Numeric vector of target recurrence periods for
#   quantile estimation.
#   fix_zeros: Logical; if TRUE, zero values are handled specially.
# Returns:
#   A list containing:
#     - eva_table: Data frame of exceedance data.
#     - statistics: Data frame of sample statistics.
#     - distr: Distribution name.
#     - method: Fitting method.
#     - params: Fitted distribution parameters.
#     - model_pred: Model predictions at observed recurrence periods.
#     - model_quant: Data frame of model quantiles at model_rperiods.
#     - target_quant: Data frame of model quantiles at target_rperiods.
#     - metrics: Data frame of goodness-of-fit metrics.
run_probmodel <- function(data,
                          method,
                          distr,
                          model_rperiods,
                          target_rperiods,
                          fix_zeros) {
  name <- paste(distr, method, sep = "_")
  model_pexc <- 1 / model_rperiods
  eva_table <- build_eva_table(data, remove_zeros = fix_zeros)
  statistics <- get_sample_statistics(data, remove_zeros = fix_zeros)
  params <- fit_probmodel(
    eva_table$data, distr = distr, method = method
  )
  model_quant <- get_model_quant(model_pexc, statistics,
                                 distr, params, fix_zeros)
  model_pred <- approx(x = model_rperiods,
                       y = model_quant,
                       xout = eva_table$rperiod)$y
  target_quant <- approx(x = model_rperiods,
                         y = model_quant,
                         xout = target_rperiods)$y
  metrics <- gofmetrics(eva_table$data, model_pred,
                        distr = distr, params = params)

  model_quant <- data.frame(model_quant, row.names = model_rperiods)
  target_quant <- data.frame(target_quant, row.names = target_rperiods)
  colnames(model_quant) <- name
  colnames(target_quant) <- name
  colnames(metrics) <- name
  list(
    eva_table = eva_table,
    statistics = statistics,
    distr = distr,
    method = method,
    params = params,
    model_pred = model_pred,
    model_quant = model_quant,
    target_quant = target_quant,
    metrics = metrics
  )
}


## Removed unused helper pad_and_bind (not referenced anywhere)