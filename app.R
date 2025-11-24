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
    ci_results = NULL
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
      h5("(Non-zero values only):"),
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
    cat("Distribution:", input$distribution, "\n")
    cat("Method:", input$method, "\n\n")
    cat("Parameters:\n")
    print(round(rv$initial_params, 3))
  })
  
  output$param_controls <- renderUI({
    req(rv$fitted_params)
    create_param_controls_ui(rv$fitted_params, input$distribution)
  })
  
  # ============= Recompute with Manual Parameters ============= #
  observeEvent(input$recompute, {
    req(rv$results)
    req(rv$fitted_params)
    
    rv$ci_results <- NULL
    
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
  output$gof_tests <- renderTable({
    req(rv$results)
    req(rv$data)
    run_gof_tests(rv$data, rv$results$distr, as.numeric(rv$results$params))
  }, rownames = FALSE, align = 'c')
  
output$quantiles_table <- renderUI({
  req(rv$results)
  qt <- format_quantiles_table(rv$results$target_quant, rv$ci_results)
  # Return periods as numeric
  periods <- rownames(qt)
  values <- as.numeric(qt[, 1])
  lower <- if (!is.null(rv$ci_results) && !is.null(qt$`Lower CI`)) as.numeric(qt$`Lower CI`) else NULL
  upper <- if (!is.null(rv$ci_results) && !is.null(qt$`Upper CI`)) as.numeric(qt$`Upper CI`) else NULL
  fmt_period <- function(x) {
    v <- as.numeric(x)
    if (is.na(v) || !is.finite(v)) return("—")
    as.character(round(v))
  }
  fmt <- function(x) {
    v <- as.numeric(x)
    if (is.na(v) || !is.finite(v)) return("—")
    format(round(v, 3), nsmall = 3)
  }
  tags$table(
    class = "table table-bordered table-hover",
    style = "text-align: center; margin: auto; max-width: 800px;",
    tags$thead(
      tags$tr(
        tags$th("Return Period", style = "text-align: center;"),
        tags$th("Return Level", style = "text-align: center;"),
        if (!is.null(lower)) tags$th("Lower CI", style = "text-align: center;"),
        if (!is.null(upper)) tags$th("Upper CI", style = "text-align: center;")
      )
    ),
    tags$tbody(
      lapply(seq_along(periods), function(i)
        tags$tr(
          tags$td(fmt_period(periods[i])),
          tags$td(fmt(values[i])),
          if (!is.null(lower)) tags$td(fmt(lower[i])),
          if (!is.null(upper)) tags$td(fmt(upper[i]))
        )
      )
    )
  )
})  # ============= Plot Outputs ============= #
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
    create_qq_plot(rv$results$eva_table, rv$results$distr, rv$results$params)
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
    create_cdf_plot(rv$results$eva_table, rv$results$distr, rv$results$params)
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
}

# ==================== Run Application ==================== #
shinyApp(ui = ui, server = server)
