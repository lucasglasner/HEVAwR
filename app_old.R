library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(openxlsx)
library(tibble)

source("fit_utils.R")
source("test_utils.R")
source("global_utils.R")

# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Extreme Value Analysis - Distribution Fitting"),
  sidebarLayout(
    # Sidebar for inputs
    sidebarPanel(
      width = 3,
      # Data upload section
      h4("1. Data Input"),
      fileInput("datafile", "Upload CSV/TXT File",
                accept = c(".csv", ".txt")),
      helpText("Upload a file with a single column of numeric values"),
      # Manual data input option
      textAreaInput("manual_data", "Or paste data (one value per line):",
                    rows = 5, placeholder = "10.5\n12.3\n15.7\n..."),
      checkboxInput("fix_zeros", "Handle zero values", value = FALSE),
      hr(),
      # Distribution selection
      h4("2. Distribution Settings"),
      selectInput("distribution", "Select Distribution:",
                  choices = c("Normal" = "norm",
                            "Lognormal" = "lognorm",
                            "Gamma" = "gamma",
                            "Pearson Type III" = "pearson3",
                            "Log-Pearson Type III" = "logpearson3",
                            "Gumbel" = "gumbel",
                            "GEV" = "gev"),
                  selected = "gev"),
      
      selectInput("method", "Fitting Method:",
                  choices = c("L-moments" = "lmme",
                            "Method of Moments" = "mme",
                            "Maximum Likelihood" = "mle"),
                  selected = "lmme"),
      
      hr(),
      
      # Return periods
      h4("3. Return Periods"),
      textInput("target_rperiods", "Target return periods (comma-separated):",
                value = "2, 5, 10, 25, 50, 100, 500"),
      
      hr(),
      
      actionButton("run_analysis", "Run Analysis", 
                   class = "btn-primary btn-block",
                   icon = icon("play"))
    ),
    
    # Main panel for outputs
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        # Data preview tab
        tabPanel("Data Preview",
                 h4("Uploaded Data"),
                 verbatimTextOutput("data_summary"),
                 DTOutput("data_table"),
                 hr(),
                 h4("Sample Statistics"),
                 tableOutput("sample_stats")
        ),
        
        # Fitting results tab
        tabPanel("Fitting Results",
                 fluidRow(
                   column(12,
                          h4("Fitted Parameters (Initial Fit)"),
                          verbatimTextOutput("initial_params"),
                          hr(),
                          h4("Manual Parameter Adjustment"),
                          helpText("Adjust parameters manually and click 'Recompute' to update results"),
                          uiOutput("param_controls"),
                          actionButton("recompute", "Recompute with Manual Parameters", 
                                      class = "btn-warning",
                                      icon = icon("calculator")),
                          hr(),
                          h4("Confidence Intervals"),
                          helpText("Compute bootstrap confidence intervals for uncertainty quantification"),
                          fluidRow(
                            column(4,
                                   numericInput("ci_level", "Confidence Level:", 
                                               value = 0.95, min = 0.5, max = 0.99, step = 0.01)
                            ),
                            column(4,
                                   numericInput("n_bootstrap", "Bootstrap Iterations:", 
                                               value = 100, min = 100, max = 10000, step = 100)
                            ),
                            column(4,
                                   checkboxInput("parallel_bootstrap", "Use Parallel", value = FALSE)
                            )
                          ),
                          actionButton("compute_ci", "Compute Confidence Intervals", 
                                      class = "btn-info",
                                      icon = icon("chart-line")),
                          hr()
                   )
                 ),
                 h4("Probability Plot (Value vs Return Period)"),
                 downloadButton("download_prob_rperiod", "Download Plot", class = "btn-sm"),
                 plotOutput("prob_plot_rperiod", height = "500px"),
                 hr(),
                 fluidRow(
                   column(6,
                          h4("Goodness-of-Fit Tests"),
                          tableOutput("gof_tests")
                   ),
                   column(6,
                          h4("Return Period Quantiles"),
                          tableOutput("quantiles_table")
                   )
                 ),
                 hr(),
                 downloadButton("download_report", "Download Excel Report", 
                               class = "btn-success")
        ),
        
        # Plots tab
        tabPanel("Plots",
                 fluidRow(
                   column(6,
                          h4("Probability Plot (Value vs Exceedance Probability)"),
                          downloadButton("download_prob_pexc", "Download Plot", class = "btn-sm"),
                          plotOutput("prob_plot", height = "400px")
                   ),
                   column(6,
                          h4("Q-Q Plot"),
                          downloadButton("download_qq", "Download Plot", class = "btn-sm"),
                          plotOutput("qq_plot", height = "400px")
                   )
                 ),
                 hr(),
                 fluidRow(
                   column(6,
                          h4("Histogram with Fitted PDF"),
                          downloadButton("download_hist", "Download Plot", class = "btn-sm"),
                          plotOutput("hist_plot", height = "400px")
                   ),
                   column(6,
                          h4("Empirical vs Fitted CDF"),
                          downloadButton("download_cdf", "Download Plot", class = "btn-sm"),
                          plotOutput("cdf_plot", height = "400px")
                   )
                 )
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive values to store data and results
  rv <- reactiveValues(
    data = NULL,
    results = NULL,
    fitted_params = NULL,
    initial_params = NULL,
    ci_results = NULL
  )
  
  # Load data from file or manual input
  observe({
    if (!is.null(input$datafile)) {
      df <- read.csv(input$datafile$datapath, header = FALSE)
      rv$data <- df[, 1]
    } else if (input$manual_data != "") {
      rv$data <- as.numeric(unlist(strsplit(input$manual_data, "\n")))
      rv$data <- rv$data[!is.na(rv$data)]
    }
  })
  
  # Data summary
  output$data_summary <- renderPrint({
    req(rv$data)
    cat("Number of observations:", length(rv$data), "\n")
    cat("Range:", paste(range(rv$data, na.rm = TRUE), collapse = " to "), "\n")
    cat("Missing values:", sum(is.na(rv$data)), "\n")
  })
  
  # Data table
  output$data_table <- renderDT({
    req(rv$data)
    datatable(data.frame(Index = 1:length(rv$data), Value = rv$data),
              options = list(pageLength = 10))
  })
  
  # Sample statistics
  output$sample_stats <- renderTable({
    req(rv$data)
    get_sample_statistics(rv$data, remove_zeros = input$fix_zeros)
  }, rownames = FALSE)
  
  # Clear CI results when distribution or method changes
  observeEvent(input$distribution, {
    rv$ci_results <- NULL
  }, ignoreInit = TRUE)
  
  observeEvent(input$method, {
    rv$ci_results <- NULL
  }, ignoreInit = TRUE)
  
  # Run main analysis
  observeEvent(input$run_analysis, {
    req(rv$data)
    
    # Clear previous CI results
    rv$ci_results <- NULL
    
    # Parse target return periods
    target_rp <- as.numeric(unlist(strsplit(input$target_rperiods, ",")))
    model_rp <- exp(seq(log(1.01), log(max(target_rp) * 2), length.out = 100))
    
    withProgress(message = 'Fitting distribution...', value = 0, {
      
      # Run the probability model
      rv$results <- run_probmodel(
        data = rv$data,
        method = input$method,
        distr = input$distribution,
        model_rperiods = model_rp,
        target_rperiods = target_rp,
        fix_zeros = input$fix_zeros
      )
      
      # Store fitted parameters
      rv$fitted_params <- rv$results$params
      rv$initial_params <- rv$results$params
    })
  })
  
  # Compute confidence intervals using bootstrap
  observeEvent(input$compute_ci, {
    req(rv$data)
    req(rv$results)
    
    withProgress(message = 'Computing confidence intervals...', value = 0, {
      # Run bootstrap
      bootstrap_results <- bootstrap_probmodel(
        x = rv$data,
        niters = input$n_bootstrap,
        distr = input$distribution,
        method = input$method,
        parallel = input$parallel_bootstrap,
        replace = TRUE,
        seed = 123,
        progress = FALSE
      )
      
      if (is.null(bootstrap_results)) {
        showNotification(
          "Bootstrap failed. Try different settings.",
          type = "error"
        )
        return()
      }
      
      setProgress(0.5, detail = "Calculating quantiles...")
      
      # Calculate confidence intervals for parameters
      alpha <- 1 - input$ci_level
      
      # Get quantiles for model predictions
      target_rp <- as.numeric(unlist(strsplit(input$target_rperiods, ",")))
      model_rp <- exp(seq(log(1.01), log(max(target_rp) * 2), length.out = 100))
      model_pexc <- 1 / model_rp
      
      # Compute quantiles for each bootstrap iteration
      n_iters <- nrow(bootstrap_results)
      model_quant_boot <- matrix(NA, nrow = length(model_rp), ncol = n_iters)
      target_quant_boot <- matrix(NA, nrow = length(target_rp), ncol = n_iters)
      
      for (i in 1:n_iters) {
        params_i <- as.numeric(bootstrap_results[i, ])
        if (!any(is.na(params_i))) {
          model_quant_boot[, i] <- get_model_quant(
            model_pexc, rv$results$statistics,
            input$distribution, params_i, input$fix_zeros
          )
          target_quant_boot[, i] <- approx(x = model_rp,
                                           y = model_quant_boot[, i],
                                           xout = target_rp)$y
        }
      }
      
      setProgress(0.8, detail = "Computing intervals...")
      
      # Compute confidence intervals
      model_quant_lower <- apply(model_quant_boot, 1, quantile, 
                                 probs = alpha/2, na.rm = TRUE)
      model_quant_upper <- apply(model_quant_boot, 1, quantile, 
                                 probs = 1 - alpha/2, na.rm = TRUE)
      
      target_quant_lower <- apply(target_quant_boot, 1, quantile, 
                                  probs = alpha/2, na.rm = TRUE)
      target_quant_upper <- apply(target_quant_boot, 1, quantile, 
                                  probs = 1 - alpha/2, na.rm = TRUE)
      
      # Store results
      rv$ci_results <- list(
        model_rp = model_rp,
        model_lower = model_quant_lower,
        model_upper = model_quant_upper,
        target_rp = target_rp,
        target_lower = target_quant_lower,
        target_upper = target_quant_upper,
        ci_level = input$ci_level
      )
      
      setProgress(1, detail = "Done!")
      showNotification(
        paste0(input$ci_level * 100, "% CI computed!"),
        type = "message"
      )
    })
  })
  
  # Display initial fitted parameters
  output$initial_params <- renderPrint({
    req(rv$initial_params)
    cat("Distribution:", input$distribution, "\n")
    cat("Method:", input$method, "\n\n")
    cat("Parameters:\n")
    print(round(rv$initial_params, 3))
  })
  
  # Dynamic parameter controls based on fitted parameters
  output$param_controls <- renderUI({
    req(rv$fitted_params)
    
    distr <- input$distribution
    n_params <- length(rv$fitted_params)
    
    param_names <- switch(distr,
      "norm" = c("Mean", "SD"),
      "lognorm" = c("Mean (log)", "SD (log)"),
      "gamma" = c("Shape", "Rate"),
      "pearson3" = c("Location", "Scale", "Shape"),
      "logpearson3" = c("Location (log)", "Scale (log)", "Shape (log)"),
      "gumbel" = c("Location", "Scale"),
      "gev" = c("Location", "Scale", "Shape")
    )
    
    fluidRow(
      lapply(1:n_params, function(i) {
        column(4,
          numericInput(paste0("manual_param_", i), 
                      param_names[i], 
                      value = round(rv$fitted_params[i], 3),
                      step = 0.001)
        )
      })
    )
  })
  
  # Recompute with manual parameters
  observeEvent(input$recompute, {
    req(rv$results)
    req(rv$fitted_params)
    
    # Clear previous CI results since parameters changed
    rv$ci_results <- NULL
    
    n_params <- length(rv$fitted_params)
    manual_params <- sapply(1:n_params, function(i) {
      input[[paste0("manual_param_", i)]]
    })
    
    # Update parameters
    rv$results$params <- manual_params
    
    # Recompute model predictions
    target_rp <- as.numeric(unlist(strsplit(input$target_rperiods, ",")))
    model_rp <- exp(seq(log(1.01), log(max(target_rp) * 2), length.out = 100))
    model_pexc <- 1 / model_rp
    
    model_quant <- get_model_quant(model_pexc, rv$results$statistics,
                                   rv$results$distr, manual_params, 
                                   input$fix_zeros)
    
    model_pred <- approx(x = model_rp,
                         y = model_quant,
                         xout = rv$results$eva_table$rperiod)$y
    
    target_quant <- approx(x = model_rp,
                           y = model_quant,
                           xout = target_rp)$y
    
    # Recompute metrics
    metrics <- gofmetrics(rv$results$eva_table$data, model_pred,
                          distr = rv$results$distr, params = manual_params)
    
    # Update results
    model_quant_df <- data.frame(model_quant, row.names = model_rp)
    target_quant_df <- data.frame(target_quant, row.names = target_rp)
    name <- paste(rv$results$distr, rv$results$method, sep = "_")
    colnames(model_quant_df) <- name
    colnames(target_quant_df) <- name
    colnames(metrics) <- name
    
    rv$results$model_quant <- model_quant_df
    rv$results$target_quant <- target_quant_df
    rv$results$metrics <- metrics
    rv$results$model_pred <- model_pred
  })
  
  # Display GOF metrics
  output$gof_metrics <- renderTable({
    req(rv$results)
    t(rv$results$metrics)
  }, rownames = TRUE)
  
  # Display GOF tests
  output$gof_tests <- renderTable({
    req(rv$results)
    req(rv$data)
    
    # Get current parameters and data
    x <- rv$data
    distr <- rv$results$distr
    # Use rv$results$params to reflect manual changes
    params <- as.numeric(rv$results$params)
    
    # Run individual GOF tests from test_utils.R
    chi2 <- chi2_gof(x, distr, params, nbins = NULL, alpha = 0.05)
    ks <- ks_gof(x, distr, params, alpha = 0.05)
    cvm <- cvm_gof(x, distr, params, alpha = 0.05)
    ad <- ad_gof(x, distr, params, alpha = 0.05)
    
    # Build results dataframe
    tests_df <- data.frame(
      Test = c(
        "Chi-Squared",
        "Kolmogorov-Smirnov",
        "Cramer-von Mises",
        "Anderson-Darling"
      ),
      Statistic = c(
        round(chi2$statistic, 4),
        round(ks$statistic, 4),
        round(cvm$statistic, 4),
        round(ad$statistic, 4)
      ),
      `P-Value` = c(
        round(chi2$pvalue, 4),
        round(ks$pvalue, 4),
        round(cvm$pvalue, 4),
        round(ad$pvalue, 4)
      ),
      Passed = c(
        ifelse(chi2$test, "✓", "✗"),
        ifelse(ks$test, "✓", "✗"),
        ifelse(cvm$test, "✓", "✗"),
        ifelse(ad$test, "✓", "✗")
      ),
      check.names = FALSE
    )
    
    tests_df
  }, rownames = FALSE)
  
  # Display quantiles with confidence intervals
  output$quantiles_table <- renderTable({
    req(rv$results)
    
    quant_df <- rv$results$target_quant
    
    # Add confidence intervals if available
    if (!is.null(rv$ci_results)) {
      ci_level_pct <- rv$ci_results$ci_level * 100
      colnames(quant_df) <- "Estimate"
      quant_df <- cbind(
        quant_df,
        `Lower CI` = rv$ci_results$target_lower,
        `Upper CI` = rv$ci_results$target_upper
      )
    }
    
    quant_df
  }, rownames = TRUE)
  
  # Probability plot (exceedance probability)
  output$prob_plot <- renderPlot({
    req(rv$results)
    
    eva_table <- rv$results$eva_table
    model_rperiods <- as.numeric(rownames(rv$results$model_quant))
    model_quant <- rv$results$model_quant[, 1]
    model_pexc <- 1 / model_rperiods
    
    ggplot() +
      geom_point(data = eva_table, 
                aes(x = pexc, y = data), 
                color = "blue", size = 3, alpha = 0.6) +
      geom_line(aes(x = model_pexc, y = model_quant), 
               color = "red", linewidth = 1.0) +
      scale_x_log10() +
      labs(x = "Exceedance Probability", 
           y = "Return Level",
           title = paste("Probability Plot -", rv$results$distr)) +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  })
  
  # Probability plot (return period) - for Fitting Results tab
  output$prob_plot_rperiod <- renderPlot({
    req(rv$results)
    
    eva_table <- rv$results$eva_table
    model_rperiods <- as.numeric(rownames(rv$results$model_quant))
    model_quant <- rv$results$model_quant[, 1]
    
    # Extract metrics
    metrics <- rv$results$metrics
    n <- round(metrics["n", 1], 0)
    r2 <- round(metrics["r2", 1], 3)
    rmse <- round(metrics["rmse", 1], 3)
    mbias <- round(metrics["mbias", 1], 3)
    aic <- round(metrics["aic", 1], 2)
    bic <- round(metrics["bic", 1], 2)
    
    # Create metrics label
    metrics_label <- paste0(
      "n = ", n, "\n",
      "R² = ", r2, "\n",
      "RMSE = ", rmse, "\n",
      "MBias = ", mbias, "\n",
      "AIC = ", aic, "\n",
      "BIC = ", bic
    )
    
    p <- ggplot() +
      geom_point(data = eva_table, 
                aes(x = rperiod, y = data), 
                color = "blue", size = 3, alpha = 0.6) +
      geom_line(aes(x = model_rperiods, y = model_quant), 
               color = "red", linewidth = 1.0) +
      scale_x_log10() +
      labs(x = "Return Period (years)", 
           y = "Return Level",
           title = paste("Probability Plot -", rv$results$distr)) +
      annotate("text", x = Inf, y = -Inf, label = metrics_label,
               hjust = 1.1, vjust = -0.1, size = 5, family = "mono",
               fontface = "plain") +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
    
    # Add confidence intervals if available
    if (!is.null(rv$ci_results)) {
      ci_data <- data.frame(
        rperiod = rv$ci_results$model_rp,
        lower = rv$ci_results$model_lower,
        upper = rv$ci_results$model_upper
      )
      p <- p +
        geom_line(data = ci_data, aes(x = rperiod, y = lower), 
                 color = "red", linetype = "dashed", linewidth = 0.8) +
        geom_line(data = ci_data, aes(x = rperiod, y = upper), 
                 color = "red", linetype = "dashed", linewidth = 0.8)
    }
    
    p
  })
  
  # Q-Q Plot
  output$qq_plot <- renderPlot({
    req(rv$results)
    
    eva_table <- rv$results$eva_table
    
    ggplot(data = eva_table, aes(sample = data)) +
      stat_qq(distribution = function(p) {
        qprobmodel(p, rv$results$distr, rv$results$params)
      }, size = 3) +
      stat_qq_line(distribution = function(p) {
        qprobmodel(p, rv$results$distr, rv$results$params)
      }, color = "red", linewidth = 1.0) +
      labs(x = "Theoretical Quantiles", 
           y = "Sample Quantiles",
           title = "Q-Q Plot") +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  })
  
  # Histogram with fitted PDF
  output$hist_plot <- renderPlot({
    req(rv$results)
    
    eva_table <- rv$results$eva_table
    x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
    pdf_vals <- dprobmodel(x_seq, rv$results$distr, rv$results$params)
    
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
           title = "Histogram with Fitted PDF") +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  })
  
  # Empirical vs Fitted CDF
  output$cdf_plot <- renderPlot({
    req(rv$results)
    
    eva_table <- rv$results$eva_table
    x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
    cdf_vals <- pprobmodel(x_seq, rv$results$distr, rv$results$params)
    
    ggplot() +
      stat_ecdf(data = eva_table, aes(x = data), 
               geom = "step", color = "blue", linewidth = 1.0) +
      geom_line(aes(x = x_seq, y = cdf_vals), color = "red", linewidth = 1.0) +
      labs(x = "Value", 
           y = "Cumulative Probability",
           title = "Empirical vs Fitted CDF") +
      theme_minimal(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14))
  })
  
  # Download Excel report
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("EVA_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$results)
      
      # Create workbook
      wb <- createWorkbook()
      
      # Add worksheets
      addWorksheet(wb, sheetName = "EVA_Table")
      addWorksheet(wb, sheetName = "Statistics")
      addWorksheet(wb, sheetName = "CalibrationParams")
      addWorksheet(wb, sheetName = "FitModel")
      addWorksheet(wb, sheetName = "FitResults")
      addWorksheet(wb, sheetName = "PerformanceMetrics")
      
      # Prepare parameter data
      param_df <- data.frame(
        distr = paste(rv$results$distr, rv$results$method, sep = "_"),
        loc = ifelse(
          length(rv$results$params) >= 1,
          rv$results$params[1], NA
        ),
        scale = ifelse(
          length(rv$results$params) >= 2,
          rv$results$params[2], NA
        ),
        shape = ifelse(
          length(rv$results$params) >= 3,
          rv$results$params[3], NA
        )
      )
      
      # Write data to worksheets
      writeData(wb, sheet = "EVA_Table", x = rv$results$eva_table)
      writeData(wb, sheet = "Statistics", x = rv$results$statistics)
      writeData(wb, sheet = "CalibrationParams", x = param_df)
      writeData(wb, sheet = "FitModel", 
                x = tibble::rownames_to_column(rv$results$model_quant, "T"))
      writeData(wb, sheet = "FitResults", 
                x = tibble::rownames_to_column(rv$results$target_quant, "T"))
      writeData(wb, sheet = "PerformanceMetrics", 
                x = tibble::rownames_to_column(rv$results$metrics, "metric"))
      
      # Save workbook
      saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )
  
  # Download handler for probability plot (return period)
  output$download_prob_rperiod <- downloadHandler(
    filename = function() {
      paste0("prob_plot_rperiod_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      eva_table <- rv$results$eva_table
      model_rperiods <- as.numeric(rownames(rv$results$model_quant))
      model_quant <- rv$results$model_quant[, 1]
      
      # Extract metrics
      metrics <- rv$results$metrics
      n <- round(metrics["n", 1], 0)
      r2 <- round(metrics["r2", 1], 3)
      rmse <- round(metrics["rmse", 1], 3)
      mbias <- round(metrics["mbias", 1], 3)
      aic <- round(metrics["aic", 1], 2)
      bic <- round(metrics["bic", 1], 2)
      
      # Create metrics label
      metrics_label <- paste0(
        "n = ", n, "\n",
        "R² = ", r2, "\n",
        "RMSE = ", rmse, "\n",
        "MBias = ", mbias, "\n",
        "AIC = ", aic, "\n",
        "BIC = ", bic
      )
      
      p <- ggplot() +
        geom_point(data = eva_table, 
                  aes(x = rperiod, y = data), 
                  color = "blue", size = 3, alpha = 0.6) +
        geom_line(aes(x = model_rperiods, y = model_quant), 
                 color = "red", linewidth = 1.0) +
        scale_x_log10() +
        labs(x = "Return Period (years)", 
             y = "Return Level",
             title = paste("Probability Plot -", rv$results$distr)) +
        annotate("text", x = Inf, y = -Inf, label = metrics_label,
                 hjust = 1.1, vjust = -0.1, size = 5, family = "mono",
                 fontface = "plain") +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
      
      # Add confidence intervals if available
      if (!is.null(rv$ci_results)) {
        ci_data <- data.frame(
          rperiod = rv$ci_results$model_rp,
          lower = rv$ci_results$model_lower,
          upper = rv$ci_results$model_upper
        )
        p <- p +
          geom_line(data = ci_data, aes(x = rperiod, y = lower), 
                   color = "red", linetype = "dashed", linewidth = 0.8) +
          geom_line(data = ci_data, aes(x = rperiod, y = upper), 
                   color = "red", linetype = "dashed", linewidth = 0.8)
      }
      
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Download handler for probability plot (exceedance probability)
  output$download_prob_pexc <- downloadHandler(
    filename = function() {
      paste0("prob_plot_pexc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      eva_table <- rv$results$eva_table
      model_rperiods <- as.numeric(rownames(rv$results$model_quant))
      model_quant <- rv$results$model_quant[, 1]
      model_pexc <- 1 / model_rperiods
      
      p <- ggplot() +
        geom_point(data = eva_table, 
                  aes(x = pexc, y = data), 
                  color = "blue", size = 3, alpha = 0.6) +
        geom_line(aes(x = model_pexc, y = model_quant), 
                 color = "red", linewidth = 1.0) +
        scale_x_log10() +
        labs(x = "Exceedance Probability", 
             y = "Return Level",
             title = paste("Probability Plot -", rv$results$distr)) +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
      
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Download handler for Q-Q plot
  output$download_qq <- downloadHandler(
    filename = function() {
      paste0("qq_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      eva_table <- rv$results$eva_table
      
      p <- ggplot(data = eva_table, aes(sample = data)) +
        stat_qq(distribution = function(p) {
          qprobmodel(p, rv$results$distr, rv$results$params)
        }, size = 3) +
        stat_qq_line(distribution = function(p) {
          qprobmodel(p, rv$results$distr, rv$results$params)
        }, color = "red", linewidth = 1.0) +
        labs(x = "Theoretical Quantiles", 
             y = "Sample Quantiles",
             title = "Q-Q Plot") +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
      
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Download handler for histogram
  output$download_hist <- downloadHandler(
    filename = function() {
      paste0("histogram_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      eva_table <- rv$results$eva_table
      x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
      pdf_vals <- dprobmodel(x_seq, rv$results$distr, rv$results$params)
      
      # Calculate number of bins using Freedman-Diaconis rule with minimum of 10
      nbins <- max(10, nclass.FD(eva_table$data))
      
      p <- ggplot() +
        geom_histogram(
          data = eva_table,
          aes(x = data, y = after_stat(density)),
          bins = nbins, fill = "lightblue",
          color = "black", alpha = 0.6,
          linewidth = 0.5
        ) +
        geom_line(
          aes(x = x_seq, y = pdf_vals),
          color = "red", linewidth = 1.0
        ) +
        labs(x = "Value", 
             y = "Density",
             title = "Histogram with Fitted PDF") +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
      
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Download handler for CDF plot
  output$download_cdf <- downloadHandler(
    filename = function() {
      paste0("cdf_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(rv$results)
      eva_table <- rv$results$eva_table
      x_seq <- seq(min(eva_table$data), max(eva_table$data), length.out = 200)
      cdf_vals <- pprobmodel(x_seq, rv$results$distr, rv$results$params)
      
      p <- ggplot() +
        stat_ecdf(data = eva_table, aes(x = data), 
                 geom = "step", color = "blue", linewidth = 1.0) +
        geom_line(aes(x = x_seq, y = cdf_vals), color = "red", linewidth = 1.0) +
        labs(x = "Value", 
             y = "Cumulative Probability",
             title = "Empirical vs Fitted CDF") +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, size = 18),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14))
      
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
