# UI Module for EVA Shiny App
# Contains all UI components organized by functionality

library(shiny)
library(shinythemes)

# ==================== UI Helper Functions ==================== #

# Sidebar panel with data inputs and controls.
create_sidebar_ui <- function() {
  sidebarPanel(
    width = 2,
    # Data name section
    h5(strong("Data Name")),
    textInput(
      "data_name",
      NULL,
      value = "Example",
      placeholder = "e.g., Station_1923812"
    ),
    p("Name for identifying this dataset",
      style = "font-size: 0.9em; color: #666; margin-top: -10px;"),
    hr(),
    # Data upload section
    h5(strong("Data Input")),
    fileInput(
      "datafile", NULL,
      accept = c(".csv", ".txt")
    ),
    p("Upload a file with a single column of numeric values", 
      style = "font-size: 0.9em; color: #666;"),
    # Manual data input option
    textAreaInput(
      "manual_data", "Or paste data (one value per line):",
      rows = 5, placeholder = "10.5\n12.3\n15.7\n..."
    ),
    hr(),
    checkboxInput(
      "handle_zeros",
      "Handle zero values",
      value = TRUE
    ),
    p("When checked, zero values are removed and probabilities are adjusted",
      style = "font-size: 0.85em; color: #666; margin-top: -10px;"),
    hr(),
    h5(strong("Return Periods")),
    textInput(
      "target_rperiods",
      "Target return periods (comma-separated):",
      value = "2, 5, 10, 20, 25, 50, 100, 150, 200"
    ),
    hr(),
    h5(strong("References")),
    tags$div(
      style = "font-size: 0.85em; color: #555;",
      tags$p(
        'Delignette-Muller, M.-L., Dutang, C., Pouillot, R., Denis, J.-B., & Siberchicot, A. (2025). fitdistrplus: Help to Fit of a Parametric Distribution to Non-Censored or Censored Data (Version 1.2-4) [R package]. CRAN.',
        style = "margin-bottom: 10px;"
      ),
      tags$p(
        "Wilks, D. S. (2011). Statistical methods in the atmospheric sciences (Vol. 100). Academic press.",
        style = "margin-bottom: 10px;"
      ),
      tags$p(
        "Hosking, J. R. M., & Wallis, J. R. (1997). Regional frequency analysis (p. 240).",
        style = "margin-bottom: 10px;"
      )
    )
  )
}

# Data preview tab (stats + table).
create_data_preview_tab <- function() {
  tabPanel(HTML("<i class='fa fa-table'></i> Data Preview"),
           h4("Sample Statistics", style = "text-align: center;"),
           uiOutput("sample_stats"),
           h4("Uploaded Data", style = "text-align: center;"),
           verbatimTextOutput("data_summary"),
           DTOutput("data_table"),
           hr()
  )
}

# Fitting results tab (params, plots, GOF, quantiles, export).
create_fitting_results_tab <- function() {
  tabPanel(HTML("<i class='fa fa-chart-line'></i> Fitting Tool"),
    fluidRow(
      column(12,
        h4("Distribution Settings", style = "text-align: center;"),
        fluidRow(
          column(4, offset = 2,
            selectInput("distribution", "Select Distribution:",
                        choices = c("Normal" = "norm",
                                  "Lognormal" = "lognorm",
                                  "Gamma" = "gamma",
                                  "Pearson Type III" = "pearson3",
                                  "Log-Pearson Type III" = "logpearson3",
                                  "Gumbel" = "gumbel",
                                  "GEV" = "gev"),
                        selected = "gev")
          ),
          column(4,
            selectInput("method", "Fitting Method:",
                        choices = c("L-moments" = "lmme",
                                  "Method of Moments" = "mme",
                                  "Maximum Likelihood" = "mle"),
                        selected = "lmme")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(12,
        div(
          style = "text-align: center; margin-bottom: 10px;",
          actionButton(
            "run_analysis", "Run Analysis",
            class = "btn-primary btn-md",
            icon = icon("play")
          ),
          downloadButton(
            "download_report",
            "Download Excel Report",
            class = "btn-success btn-md",
            icon = icon("file-excel")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(4,
        tags$div(
          style = "height: calc(100vh - 250px); overflow-y: auto; padding-right: 10px;",
          h4("Fitted Parameters (Initial Fit)", style = "text-align: center;"),
          verbatimTextOutput("initial_params"),
          hr(),
          h4("Goodness-of-Fit Tests", style = "text-align: center;"),
          verbatimTextOutput("gof_tests"),
          hr(),
          h4("Manual Parameter Adjustment", style = "text-align: center;"),
          p("Adjust parameters manually and click 'Recompute' to update results",
            style = "font-size: 0.9em; color: #666; text-align: center;"),
          uiOutput("param_controls"),
          div(
            style = "text-align: center;",
            actionButton(
              "recompute",
              "Recompute with Manual Parameters",
              class = "btn-warning",
              icon = icon("calculator")
            )
          ),
          hr(),
          h4("Confidence Intervals", style = "text-align: center;"),
          p("Compute bootstrap confidence intervals for uncertainty quantification",
            style = "font-size: 0.9em; color: #666; text-align: center;"),
          fluidRow(
            column(
              6,
              numericInput(
                "ci_level", "Confidence Level:",
                value = 0.95, min = 0.5,
                max = 0.99, step = 0.01
              )
            ),
            column(
              6,
              numericInput(
                "n_bootstrap", "Bootstrap Iterations:",
                value = 100, min = 100,
                max = 10000, step = 100
              )
            )
          ),
          div(
            style = "text-align: center;",
            actionButton(
              "compute_ci",
              "Compute Confidence Intervals",
              class = "btn-info",
              icon = icon("chart-line")
            )
          ),
          hr()
        )
      ),
      column(8,
        uiOutput("prob_plot_title"),
        fluidRow(
          column(3,
            downloadButton(
              "download_prob_rperiod",
              "Download Plot",
              class = "btn-sm"
            )
          ),
          column(3,
            numericInput(
              "ylim_min",
              "Y-axis Min:",
              value = NA,
              step = 1
            )
          ),
          column(3,
            numericInput(
              "ylim_max",
              "Y-axis Max:",
              value = NA,
              step = 1
            )
          ),
          column(3,
            actionButton(
              "reset_ylim",
              "Auto Y-axis",
              class = "btn-secondary btn-sm",
              style = "margin-top: 25px;"
            )
          )
        ),
        plotOutput("prob_plot_rperiod", height = "500px"),
        hr()
      )
    ),
    fluidRow(
      column(12,
        h4("Return Period Quantiles", style = "text-align: center;"),
        uiOutput("quantiles_table")
      )
    ),
    br(),
    hr()
  )
}

# Diagnostic plots tab (probability, QQ, histogram, CDF).
create_plots_tab <- function() {
  tabPanel(HTML("<i class='fa fa-chart-area'></i> Auxiliary Plots"),
    fluidRow(
      column(12,
        h4("Distribution Settings", style = "text-align: center;"),
        fluidRow(
          column(4, offset = 2,
            selectInput("aux_distribution", "Select Distribution:",
                        choices = c("Normal" = "norm",
                                  "Lognormal" = "lognorm",
                                  "Gamma" = "gamma",
                                  "Pearson Type III" = "pearson3",
                                  "Log-Pearson Type III" = "logpearson3",
                                  "Gumbel" = "gumbel",
                                  "GEV" = "gev"),
                        selected = "gev")
          ),
          column(4,
            selectInput("aux_method", "Fitting Method:",
                        choices = c("L-moments" = "lmme",
                                  "Method of Moments" = "mme",
                                  "Maximum Likelihood" = "mle"),
                        selected = "lmme")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(12,
        div(
          style = "text-align: center; margin-bottom: 10px;",
          actionButton(
            "aux_run_analysis", "Run Analysis",
            class = "btn-primary btn-md",
            icon = icon("play")
          ),
          downloadButton(
            "aux_download_report",
            "Download Excel Report",
            class = "btn-success btn-md",
            icon = icon("file-excel")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(
        6,
        h4("Probability Plot (Value vs Exceedance Probability)"),
        downloadButton(
          "download_prob_pexc",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("prob_plot", height = "400px")
      ),
      column(
        6,
        h4("Q-Q Plot"),
        downloadButton(
          "download_qq",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("qq_plot", height = "400px")
      )
    ),
    hr(),
    fluidRow(
      column(
        6,
        h4("Histogram with Fitted PDF"),
        downloadButton(
          "download_hist",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("hist_plot", height = "400px")
      ),
      column(
        6,
        h4("Empirical vs Fitted CDF"),
        downloadButton(
          "download_cdf",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("cdf_plot", height = "400px")
      )
    ),
    br(),
    hr()
  )
}

# Method comparison tab (multi-method vs one distribution).
create_method_comparison_tab <- function() {
  tabPanel(HTML("<i class='fa fa-balance-scale'></i> Method Comparison"),
    h3("Multi-Method Comparison", style = "text-align: center;"),
    fluidRow(
      column(12,
        h4("Select Distribution and Methods", style = "text-align: center;"),
        fluidRow(
          column(4, offset = 2,
            selectInput("method_comparison_distribution",
                       "Select Distribution:",
                       choices = c(
                         "Normal" = "norm",
                         "Lognormal" = "lognorm",
                         "Gamma" = "gamma",
                         "Pearson Type III" = "pearson3",
                         "Log-Pearson Type III" = "logpearson3",
                         "Gumbel" = "gumbel",
                         "GEV" = "gev"
                       ),
                       selected = "gev")
          ),
          column(4,
            h5("Fitting Methods:"),
            checkboxGroupInput(
              "compare_methods",
              NULL,
              choices = c(
                "L-moments" = "lmme",
                "Method of Moments" = "mme",
                "Maximum Likelihood" = "mle"
              ),
              selected = NULL,
              inline = FALSE
            )
          )
        ),
        tags$div(
          style = "text-align: center;",
          actionButton(
            "run_method_comparison",
            "Run Comparison",
            class = "btn-primary",
            icon = icon("chart-line")
          ),
          downloadButton(
            "download_method_comparison_report",
            "Download Excel Report",
            class = "btn-success",
            icon = icon("file-excel")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(3,
        h4("Method Parameters", style = "text-align: center;"),
        tags$div(
          style = "max-height: 450px; overflow-y: auto;",
          uiOutput("method_comparison_params_ui")
        )
      ),
      column(9,
        h4("Comparison Probability Plot", style = "text-align: center;"),
        downloadButton(
          "download_method_comparison_plot",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("method_comparison_plot", height = "450px")
      )
    ),
    hr(),
    fluidRow(
      column(12,
        h4("Return Period Quantiles", style = "text-align: center;"),
        uiOutput("method_comparison_quantiles_table")
      )
    ),
    br(),
    hr()
  )
}

# Model comparison tab (multi-distribution same method).
create_model_comparison_tab <- function() {
  tabPanel(HTML("<i class='fa fa-layer-group'></i> Model Comparison"),
    h3("Multi-Distribution Comparison", style = "text-align: center;"),
    fluidRow(
      column(12,
        h4("Select Distributions to Compare", style = "text-align: center;"),
        fluidRow(
          column(6, offset = 3,
            checkboxGroupInput(
              "compare_distributions",
              NULL,
              choices = c(
                "Normal" = "norm",
                "Lognormal" = "lognorm",
                "Gamma" = "gamma",
                "Pearson Type III" = "pearson3",
                "Log-Pearson Type III" = "logpearson3",
                "Gumbel" = "gumbel",
                "GEV" = "gev"
              ),
              selected = NULL,
              inline = TRUE
            )
          )
        ),
        fluidRow(
          column(4, offset = 4,
            selectInput("comparison_method", "Fitting Method:",
                       choices = c("L-moments" = "lmme",
                                 "Method of Moments" = "mme",
                                 "Maximum Likelihood" = "mle"),
                       selected = "lmme")
          )
        ),
        tags$div(
          style = "text-align: center;",
          actionButton(
            "run_comparison",
            "Run Comparison",
            class = "btn-primary",
            icon = icon("chart-line")
          ),
          downloadButton(
            "download_model_comparison_report",
            "Download Excel Report",
            class = "btn-success",
            icon = icon("file-excel")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(3,
        h4("Distribution Parameters", style = "text-align: center;"),
        tags$div(
          style = "max-height: 450px; overflow-y: auto;",
          uiOutput("comparison_params_ui")
        )
      ),
      column(9,
        h4("Comparison Probability Plot", style = "text-align: center;"),
        downloadButton(
          "download_comparison_plot",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("comparison_plot", height = "450px")
      )
    ),
    hr(),
    fluidRow(
      column(12,
        h4("Return Period Quantiles", style = "text-align: center;"),
        uiOutput("comparison_quantiles_table")
      )
    ),
    br(),
    hr()
  )
}

# Main application UI layout.
create_ui <- function() {
  fluidPage(
    theme = shinytheme("flatly"),
    titlePanel(
      title = div("Extreme Value Analysis - Distribution Fitting Tool",
          style = "text-align: center;"),
      windowTitle = "Extreme Value Analysis - Distribution Fitting Tool"
    ),
    sidebarLayout(
      create_sidebar_ui(),
      mainPanel(
        width = 10,
        tabsetPanel(
          id = "main_tabs",
          create_data_preview_tab(),
          create_fitting_results_tab(),
          create_plots_tab(),
          create_model_comparison_tab(),
          create_method_comparison_tab()
        )
      )
    )
  )
}
