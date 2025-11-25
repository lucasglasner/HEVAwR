# EVA Shiny App - Main Application File
# Extreme Value Analysis with Distribution Fitting
# Modular version with organized code structure

# ==================== Load Required Libraries ==================== #
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(openxlsx)
library(tibble)

# ==================== Source Application Modules ==================== #
source("fit_utils.R")
source("test_utils.R")
source("global_utils.R")
source("app_ui.R")
source("app_server_functions.R")
source("app_plot_functions.R")

# ==================== Define UI ==================== #
ui <- create_ui()

# ==================== Define Server ==================== #
server <- function(input, output, session) {
  
  # ============= Reactive Values ============= #
  rv <- reactiveValues(
    data = NULL,
    results = NULL,
    fitted_params = NULL,
    initial_params = NULL,
    ci_results = NULL,
    comparison_results = list(),
    method_comparison_results = list()
  )
  
  # ============= Data Loading ============= #
  observe({
    if (!is.null(input$datafile)) {
      rv$data <- load_data_from_file(input$datafile$datapath)
    } else if (input$manual_data != "") {
      rv$data <- load_data_from_text(input$manual_data)
    }
  })
  
  # ============= Data Preview Outputs ============= #
  output$data_summary <- renderPrint({
    req(rv$data)
    cat("Number of observations:", length(rv$data), "\n")
    cat("Non-zero values:", sum(rv$data != 0, na.rm = TRUE), "\n")
    cat("Missing values:", sum(is.na(rv$data)), "\n")
  })
  
  output$data_table <- renderDT({
    req(rv$data)
    n <- length(rv$data)
    sorted_data <- sort(rv$data, decreasing = TRUE)
    empirical_prob <- (1:n) / (n + 1)
    return_period <- 1 / empirical_prob
    
    df <- data.frame(
      Index = 1:n,
      Value = sorted_data,
      `P(exc)` = round(empirical_prob, 3),
      `T(years)` = round(return_period, 1),
      check.names = FALSE
    )
    
    datatable(
      df,
      options = list(
        pageLength = 5,
        columnDefs = list(
          list(className = 'dt-center', targets = '_all')
        )
      ),
      rownames = FALSE
    )
  })
  
  output$sample_stats <- renderUI({
    req(rv$data)
    stats <- get_sample_statistics(rv$data, remove_zeros = TRUE)
    fmt <- function(x) {
      v <- as.numeric(x)
      if (is.na(v) || !is.finite(v)) return("—")
      format(round(v, 3), nsmall = 3)
    }
    tagList(
      h5("(Non-zero values only)", 
         style = "text-align: center;"),
      tags$table(
        class = "table table-bordered table-hover",
        style = "text-align: center; margin: auto; max-width: 800px;",
        tags$thead(
          tags$tr(
            tags$th("Mean", style = "text-align: center;"),
            tags$th("Std Dev", style = "text-align: center;"),
            tags$th("Skewness", style = "text-align: center;"),
            tags$th("Kurtosis", style = "text-align: center;"),
            tags$th("Minimum", style = "text-align: center;"),
            tags$th("Maximum", style = "text-align: center;"),
            tags$th("Non-Zero Prob", style = "text-align: center;")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td(fmt(stats$mean)),
            tags$td(fmt(stats$sd)),
            tags$td(fmt(stats$skew)),
            tags$td(fmt(stats$kurt)),
            tags$td(fmt(stats$min)),
            tags$td(fmt(stats$max)),
            tags$td(fmt(stats$pnonzero))
          )
        )
      )
    )
  })
  
  # ============= Clear CI on Changes ============= #
  observeEvent(input$distribution, {
    rv$ci_results <- NULL
  }, ignoreInit = TRUE)
  
  observeEvent(input$method, {
    rv$ci_results <- NULL
  }, ignoreInit = TRUE)
  
  # ============= Main Analysis ============= #
  observeEvent(input$run_analysis, {
    req(rv$data)
    
    rv$ci_results <- NULL
    
    withProgress(message = 'Fitting distribution...', value = 0, {
      rv$results <- run_eva_analysis(
        data = rv$data,
        method = input$method,
        distr = input$distribution,
        target_rperiods = input$target_rperiods,
        fix_zeros = TRUE
      )
      
      rv$fitted_params <- rv$results$params
      rv$initial_params <- rv$results$params
    })
  })
  
  # ============= Parameter Controls ============= #
  output$initial_params <- renderPrint({
    req(rv$initial_params)
    req(rv$results)
    cat("Distribution:", input$distribution, "\n")
    cat("Method:", input$method, "\n\n")
    cat("Parameters:\n")
    print(round(rv$initial_params, 3))
  })
  
  output$gof_tests <- renderPrint({
    req(rv$results)
    gof <- run_gof_tests(rv$data, rv$results$distr, 
                         as.numeric(rv$results$params))
    for (i in 1:nrow(gof)) {
      cat(sprintf("%s: stat=%.4f, p=%.4f %s\n",
                  gof$Test[i], gof$Statistic[i], 
                  gof$`P-Value`[i], gof$Passed[i]))
    }
  })
  
  output$param_controls <- renderUI({
    req(rv$fitted_params)
    create_param_controls_ui(rv$fitted_params, input$distribution)
  })
  
  # ============= Recompute with Manual Parameters ============= #
  observeEvent(input$recompute, {
    req(rv$results)
    req(rv$fitted_params)
    
    # Keep CI results as reference - don't clear them
    # CIs from fitted parameters serve as uncertainty bounds guide
    
    n_params <- length(rv$fitted_params)
    manual_params <- sapply(1:n_params, function(i) {
      input[[paste0("manual_param_", i)]]
    })
    
    rv$results <- recompute_with_manual_params(
      results = rv$results,
      manual_params = manual_params,
      target_rperiods = input$target_rperiods,
      fix_zeros = TRUE
    )
  })
  
  # ============= Confidence Intervals ============= #
  observeEvent(input$compute_ci, {
    req(rv$data)
    req(rv$results)
    
    withProgress(message = 'Computing confidence intervals...', value = 0, {
      setProgress(0.3, detail = "Running bootstrap...")
      
      ci_results <- compute_bootstrap_ci(
        data = rv$data,
        n_bootstrap = input$n_bootstrap,
        distr = input$distribution,
        method = input$method,
        ci_level = input$ci_level,
        target_rperiods = input$target_rperiods,
        statistics = rv$results$statistics,
        fix_zeros = TRUE,
        parallel = input$parallel_bootstrap
      )
      
      if (is.null(ci_results)) {
        showNotification(
          "Bootstrap failed. Try different settings.",
          type = "error"
        )
        return()
      }
      
      rv$ci_results <- ci_results
      setProgress(1, detail = "Done!")
      showNotification(
        paste0(input$ci_level * 100, "% CI computed!"),
        type = "message"
      )
    })
  })
  
  # ============= GOF Tests and Tables ============= #
  output$quantiles_table <- renderUI({
  req(rv$results)
  qt <- format_quantiles_table(rv$results$target_quant, rv$ci_results)
  # Return periods as columns
  periods <- rownames(qt)
  values <- as.numeric(qt[, 1])
  lower <- if (!is.null(rv$ci_results) && !is.null(qt$`Lower CI`)) as.numeric(qt$`Lower CI`) else NULL
  upper <- if (!is.null(rv$ci_results) && !is.null(qt$`Upper CI`)) as.numeric(qt$`Upper CI`) else NULL
  
  fmt <- function(x) {
    v <- as.numeric(x)
    if (is.na(v) || !is.finite(v)) return("—")
    format(round(v, 3), nsmall = 3)
  }
  
  # Build table header with return periods as columns
  header_row <- tags$tr(
    tags$th("", style = "text-align: center;"),
    lapply(periods, function(p) 
      tags$th(as.character(round(as.numeric(p))), 
              style = "text-align: center;"))
  )
  
  # Build return level row
  level_row <- tags$tr(
    tags$th("Return Level", style = "text-align: center;"),
    lapply(values, function(v) tags$td(fmt(v)))
  )
  
  # Build CI rows if available
  lower_row <- if (!is.null(lower)) {
    tags$tr(
      tags$th("Lower CI", style = "text-align: center;"),
      lapply(lower, function(v) tags$td(fmt(v)))
    )
  } else NULL
  
  upper_row <- if (!is.null(upper)) {
    tags$tr(
      tags$th("Upper CI", style = "text-align: center;"),
      lapply(upper, function(v) tags$td(fmt(v)))
    )
  } else NULL
  
  tags$table(
    class = "table table-bordered table-hover",
    style = "text-align: center; margin: auto;",
    tags$thead(header_row),
    tags$tbody(
      level_row,
      lower_row,
      upper_row
    )
  )
})  # ============= Plot Outputs ============= #
  output$prob_plot_title <- renderUI({
    req(rv$results)
    h4(paste("Probability Plot -", toupper(rv$results$distr)), 
       style = "text-align: center;")
  })
  
  output$prob_plot_rperiod <- renderPlot({
    req(rv$results)
    create_prob_plot_rperiod(
      eva_table = rv$results$eva_table,
      model_rperiods = as.numeric(rownames(rv$results$model_quant)),
      model_quant = rv$results$model_quant[, 1],
      distr = rv$results$distr,
      metrics = rv$results$metrics,
      ci_results = rv$ci_results
    )
  })
  
  output$prob_plot <- renderPlot({
    req(rv$results)
    model_rperiods <- as.numeric(rownames(rv$results$model_quant))
    create_prob_plot_pexc(
      eva_table = rv$results$eva_table,
      model_pexc = 1 / model_rperiods,
      model_quant = rv$results$model_quant[, 1],
      distr = rv$results$distr
    )
  })
  
  output$qq_plot <- renderPlot({
    req(rv$results)
    create_qq_plot(rv$results$eva_table, rv$results$distr, 
                   rv$results$params)
  })
  
  output$hist_plot <- renderPlot({
    req(rv$results)
    create_histogram_plot(
      rv$results$eva_table,
      rv$results$distr,
      rv$results$params
    )
  })
  
  output$cdf_plot <- renderPlot({
    req(rv$results)
    create_cdf_plot(rv$results$eva_table, rv$results$distr, 
                    rv$results$params)
  })
  
  # ============= Download Handlers ============= #
  
  # Excel Report
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("EVA_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$results)
      create_excel_report(rv$results, file, rv$ci_results, rv$initial_params)
    }
  )
  
  # Probability Plot (Return Period)
  output$download_prob_rperiod <- downloadHandler(
    filename = function() {
      paste0("prob_plot_rperiod_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      p <- create_prob_plot_rperiod(
        eva_table = rv$results$eva_table,
        model_rperiods = as.numeric(rownames(rv$results$model_quant)),
        model_quant = rv$results$model_quant[, 1],
        distr = rv$results$distr,
        metrics = rv$results$metrics,
        ci_results = rv$ci_results
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Probability Plot (Exceedance Probability)
  output$download_prob_pexc <- downloadHandler(
    filename = function() {
      paste0("prob_plot_pexc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      model_rperiods <- as.numeric(rownames(rv$results$model_quant))
      p <- create_prob_plot_pexc(
        eva_table = rv$results$eva_table,
        model_pexc = 1 / model_rperiods,
        model_quant = rv$results$model_quant[, 1],
        distr = rv$results$distr
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Q-Q Plot
  output$download_qq <- downloadHandler(
    filename = function() {
      paste0("qq_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      p <- create_qq_plot(
        rv$results$eva_table,
        rv$results$distr,
        rv$results$params
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Histogram
  output$download_hist <- downloadHandler(
    filename = function() {
      paste0("histogram_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      p <- create_histogram_plot(
        rv$results$eva_table,
        rv$results$distr,
        rv$results$params
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # CDF Plot
  output$download_cdf <- downloadHandler(
    filename = function() {
      paste0("cdf_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      p <- create_cdf_plot(
        rv$results$eva_table,
        rv$results$distr,
        rv$results$params
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # ============= Model Comparison ============= #
  
  # Run comparison analysis for selected distributions
  observeEvent(input$run_comparison, {
    req(rv$data)
    req(input$compare_distributions)
    
    # Clear previous results
    rv$comparison_results <- list()
    
    # Get the common fitting method
    method <- input$comparison_method
    
    # Run analysis for each selected distribution
    for (distr in input$compare_distributions) {
      result <- run_eva_analysis(
        data = rv$data,
        method = method,
        distr = distr,
        target_rperiods = input$target_rperiods,
        fix_zeros = TRUE
      )
      rv$comparison_results[[distr]] <- result
    }
  })
  
  # Generate parameter UI for all compared distributions
  output$comparison_params_ui <- renderUI({
    req(rv$comparison_results)
    req(length(rv$comparison_results) > 0)
    
    # Create a vertical list of distribution parameter sections
    param_sections <- lapply(names(rv$comparison_results), function(distr) {
      result <- rv$comparison_results[[distr]]
      param_names <- get_param_names(distr)
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$div(
          style = "padding: 8px; border: 1px solid #ddd; 
                   border-radius: 5px; background-color: #f9f9f9;",
          h6(toupper(distr), style = "text-align: center; 
                                      font-weight: bold; 
                                      margin-top: 0px;
                                      margin-bottom: 8px;
                                      font-size: 13px;"),
          lapply(seq_along(param_names), function(i) {
            param_id <- paste0("comp_", distr, "_param", i)
            # First two parameters (location, scale) use 2 decimals
            # Other parameters (shape, etc.) use 4 decimals
            decimals <- if (i <= 2) 2 else 4
            step_size <- if (i <= 2) 0.01 else 0.001
            
            tags$div(
              style = "margin-bottom: 5px;",
              numericInput(
                param_id,
                paste0(param_names[i], ":"),
                value = round(as.numeric(result$params[i]), decimals),
                step = step_size
              )
            )
          }),
          actionButton(
            paste0("update_comp_", distr),
            "Update",
            class = "btn-sm btn-info btn-block",
            icon = icon("sync"),
            style = "margin-top: 5px; padding: 4px 8px; font-size: 12px;"
          )
        )
      )
    })
    
    tagList(param_sections)
  })
  
  # Handle parameter updates for each distribution
  observe({
    req(rv$comparison_results)
    
    lapply(names(rv$comparison_results), function(distr) {
      observeEvent(input[[paste0("update_comp_", distr)]], {
        param_names <- get_param_names(distr)
        new_params <- sapply(seq_along(param_names), function(i) {
          input[[paste0("comp_", distr, "_param", i)]]
        })
        
        # Recompute with manual parameters (pass full results list)
        updated_result <- recompute_with_manual_params(
          results = rv$comparison_results[[distr]],
          manual_params = new_params,
          target_rperiods = input$target_rperiods,
          fix_zeros = TRUE
        )
        
        rv$comparison_results[[distr]] <- updated_result
      })
    })
  })
  
  # Comparison plot
  output$comparison_plot <- renderPlot({
    req(rv$comparison_results)
    req(length(rv$comparison_results) > 0)
    
    create_comparison_plot(rv$comparison_results)
  })
  
  # Download comparison plot
  output$download_comparison_plot <- downloadHandler(
    filename = function() {
      paste0(
        "comparison_plot_", 
        format(Sys.time(), "%Y%m%d_%H%M%S"), 
        ".png"
      )
    },
    content = function(file) {
      req(rv$comparison_results)
      p <- create_comparison_plot(rv$comparison_results)
      ggsave(file, plot = p, width = 12, height = 7, dpi = 300)
    }
  )
  
  # Download model comparison Excel report
  output$download_model_comparison_report <- downloadHandler(
    filename = function() {
      paste0("Model_Comparison_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$comparison_results)
      req(length(rv$comparison_results) > 0)
      
      # Create workbook
      wb <- createWorkbook()
      
      # Add worksheets (same structure as main fitting tool)
      addWorksheet(wb, "EVA_Table")
      addWorksheet(wb, "Statistics")
      addWorksheet(wb, "CalibrationParams")
      addWorksheet(wb, "FitModel")
      addWorksheet(wb, "FitResults")
      addWorksheet(wb, "PerformanceMetrics")
      
      # Get first result for common data
      first_result <- rv$comparison_results[[1]]
      
      # 1. EVA_Table - same for all distributions
      writeData(wb, "EVA_Table", first_result$eva_table)
      
      # 2. Statistics - same for all distributions
      writeData(wb, "Statistics", first_result$statistics)
      
      # 3. CalibrationParams - all distributions as rows
      param_rows <- lapply(names(rv$comparison_results), function(distr) {
        result <- rv$comparison_results[[distr]]
        data.frame(
          distr = paste(distr, result$method, sep = "_"),
          loc = ifelse(length(result$params) >= 1, result$params[1], NA),
          scale = ifelse(length(result$params) >= 2, result$params[2], NA),
          shape = ifelse(length(result$params) >= 3, result$params[3], NA),
          stringsAsFactors = FALSE
        )
      })
      param_df <- do.call(rbind, param_rows)
      writeData(wb, "CalibrationParams", param_df)
      
      # 4. FitModel - return periods as rows, distributions as columns
      model_rp <- rownames(first_result$model_quant)
      fitmodel_df <- data.frame(T = model_rp, stringsAsFactors = FALSE)
      for (distr in names(rv$comparison_results)) {
        result <- rv$comparison_results[[distr]]
        fitmodel_df[[paste(distr, result$method, sep = "_")]] <- result$model_quant[, 1]
      }
      writeData(wb, "FitModel", fitmodel_df)
      
      # 5. FitResults - target return periods as rows, distributions as columns
      target_rp <- rownames(first_result$target_quant)
      fitresults_df <- data.frame(T = target_rp, stringsAsFactors = FALSE)
      for (distr in names(rv$comparison_results)) {
        result <- rv$comparison_results[[distr]]
        fitresults_df[[paste(distr, result$method, sep = "_")]] <- result$target_quant[, 1]
      }
      writeData(wb, "FitResults", fitresults_df)
      
      # 6. PerformanceMetrics - metrics as rows, distributions as columns
      metric_names <- rownames(first_result$metrics)
      metrics_df <- data.frame(metric = metric_names, stringsAsFactors = FALSE)
      for (distr in names(rv$comparison_results)) {
        result <- rv$comparison_results[[distr]]
        metrics_df[[paste(distr, result$method, sep = "_")]] <- result$metrics[, 1]
      }
      writeData(wb, "PerformanceMetrics", metrics_df)
      
      saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )
  
  # ============= Method Comparison ============= #
  
  # Run comparison analysis for selected methods
  observeEvent(input$run_method_comparison, {
    req(rv$data)
    req(input$compare_methods)
    
    # Clear previous results
    rv$method_comparison_results <- list()
    
    # Get the common distribution
    distr <- input$method_comparison_distribution
    
    # Run analysis for each selected method
    for (method in input$compare_methods) {
      result <- run_eva_analysis(
        data = rv$data,
        method = method,
        distr = distr,
        target_rperiods = input$target_rperiods,
        fix_zeros = TRUE
      )
      # Store with unique key (method name)
      rv$method_comparison_results[[method]] <- result
    }
  })
  
  # Generate parameter UI for all compared methods
  output$method_comparison_params_ui <- renderUI({
    req(rv$method_comparison_results)
    req(length(rv$method_comparison_results) > 0)
    
    distr <- input$method_comparison_distribution
    param_names <- get_param_names(distr)
    
    # Create a vertical list of method parameter sections
    param_sections <- lapply(names(rv$method_comparison_results), function(method) {
      result <- rv$method_comparison_results[[method]]
      
      method_label <- switch(method,
        "lmme" = "L-moments",
        "mme" = "Method of Moments",
        "mle" = "Maximum Likelihood"
      )
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$div(
          style = "padding: 8px; border: 1px solid #ddd; 
                   border-radius: 5px; background-color: #f9f9f9;",
          h6(method_label, style = "text-align: center; 
                                      font-weight: bold; 
                                      margin-top: 0px;
                                      margin-bottom: 8px;
                                      font-size: 13px;"),
          lapply(seq_along(param_names), function(i) {
            param_id <- paste0("mcomp_", method, "_param", i)
            # First two parameters (location, scale) use 2 decimals
            # Other parameters (shape, etc.) use 4 decimals
            decimals <- if (i <= 2) 2 else 4
            step_size <- if (i <= 2) 0.01 else 0.001
            
            tags$div(
              style = "margin-bottom: 5px;",
              numericInput(
                param_id,
                paste0(param_names[i], ":"),
                value = round(as.numeric(result$params[i]), decimals),
                step = step_size
              )
            )
          }),
          actionButton(
            paste0("update_mcomp_", method),
            "Update",
            class = "btn-sm btn-info btn-block",
            icon = icon("sync"),
            style = "margin-top: 5px; padding: 4px 8px; font-size: 12px;"
          )
        )
      )
    })
    
    tagList(param_sections)
  })
  
  # Handle parameter updates for each method
  observe({
    req(rv$method_comparison_results)
    
    distr <- input$method_comparison_distribution
    param_names <- get_param_names(distr)
    
    lapply(names(rv$method_comparison_results), function(method) {
      observeEvent(input[[paste0("update_mcomp_", method)]], {
        new_params <- sapply(seq_along(param_names), function(i) {
          input[[paste0("mcomp_", method, "_param", i)]]
        })
        
        # Recompute with manual parameters
        updated_result <- recompute_with_manual_params(
          results = rv$method_comparison_results[[method]],
          manual_params = new_params,
          target_rperiods = input$target_rperiods,
          fix_zeros = TRUE
        )
        
        rv$method_comparison_results[[method]] <- updated_result
      })
    })
  })
  
  # Method comparison plot
  output$method_comparison_plot <- renderPlot({
    req(rv$method_comparison_results)
    req(length(rv$method_comparison_results) > 0)
    
    create_method_comparison_plot(rv$method_comparison_results)
  })
  
  # Download method comparison plot
  output$download_method_comparison_plot <- downloadHandler(
    filename = function() {
      paste0(
        "method_comparison_plot_", 
        format(Sys.time(), "%Y%m%d_%H%M%S"), 
        ".png"
      )
    },
    content = function(file) {
      req(rv$method_comparison_results)
      p <- create_method_comparison_plot(rv$method_comparison_results)
      ggsave(file, plot = p, width = 12, height = 7, dpi = 300)
    }
  )
  
  # Download method comparison Excel report
  output$download_method_comparison_report <- downloadHandler(
    filename = function() {
      paste0("Method_Comparison_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$method_comparison_results)
      req(length(rv$method_comparison_results) > 0)
      
      # Create workbook
      wb <- createWorkbook()
      
      # Add worksheets (same structure as main fitting tool)
      addWorksheet(wb, "EVA_Table")
      addWorksheet(wb, "Statistics")
      addWorksheet(wb, "CalibrationParams")
      addWorksheet(wb, "FitModel")
      addWorksheet(wb, "FitResults")
      addWorksheet(wb, "PerformanceMetrics")
      
      # Get first result for common data
      first_result <- rv$method_comparison_results[[1]]
      distr <- input$method_comparison_distribution
      
      # 1. EVA_Table - same for all methods
      writeData(wb, "EVA_Table", first_result$eva_table)
      
      # 2. Statistics - same for all methods
      writeData(wb, "Statistics", first_result$statistics)
      
      # 3. CalibrationParams - all methods as rows
      param_rows <- lapply(names(rv$method_comparison_results), function(method) {
        result <- rv$method_comparison_results[[method]]
        data.frame(
          distr = paste(distr, method, sep = "_"),
          loc = ifelse(length(result$params) >= 1, result$params[1], NA),
          scale = ifelse(length(result$params) >= 2, result$params[2], NA),
          shape = ifelse(length(result$params) >= 3, result$params[3], NA),
          stringsAsFactors = FALSE
        )
      })
      param_df <- do.call(rbind, param_rows)
      writeData(wb, "CalibrationParams", param_df)
      
      # 4. FitModel - return periods as rows, methods as columns
      model_rp <- rownames(first_result$model_quant)
      fitmodel_df <- data.frame(T = model_rp, stringsAsFactors = FALSE)
      for (method in names(rv$method_comparison_results)) {
        result <- rv$method_comparison_results[[method]]
        fitmodel_df[[paste(distr, method, sep = "_")]] <- result$model_quant[, 1]
      }
      writeData(wb, "FitModel", fitmodel_df)
      
      # 5. FitResults - target return periods as rows, methods as columns
      target_rp <- rownames(first_result$target_quant)
      fitresults_df <- data.frame(T = target_rp, stringsAsFactors = FALSE)
      for (method in names(rv$method_comparison_results)) {
        result <- rv$method_comparison_results[[method]]
        fitresults_df[[paste(distr, method, sep = "_")]] <- result$target_quant[, 1]
      }
      writeData(wb, "FitResults", fitresults_df)
      
      # 6. PerformanceMetrics - metrics as rows, methods as columns
      metric_names <- rownames(first_result$metrics)
      metrics_df <- data.frame(metric = metric_names, stringsAsFactors = FALSE)
      for (method in names(rv$method_comparison_results)) {
        result <- rv$method_comparison_results[[method]]
        metrics_df[[paste(distr, method, sep = "_")]] <- result$metrics[, 1]
      }
      writeData(wb, "PerformanceMetrics", metrics_df)
      
      saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )
}

# ==================== Run Application ==================== #
shinyApp(ui = ui, server = server)
