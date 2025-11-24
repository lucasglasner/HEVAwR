# UI Module for EVA Shiny App
# Contains all UI components organized by functionality

library(shiny)
library(shinythemes)

# ==================== UI Helper Functions ==================== #

# Create the sidebar panel UI with data input and analysis controls.
# Returns:
#   A Shiny sidebarPanel object containing file upload inputs,
#   distribution selection, fitting method selection, and analysis
#   trigger button.
create_sidebar_ui <- function() {
  sidebarPanel(
    width = 3,
    # Data upload section
    h4("1. Data Input"),
    fileInput(
      "datafile", "Upload CSV/TXT File",
      accept = c(".csv", ".txt")
    ),
    helpText(
      "Upload a file with a single column of numeric values"
    ),
    # Manual data input option
    textAreaInput(
      "manual_data", "Or paste data (one value per line):",
      rows = 5, placeholder = "10.5\n12.3\n15.7\n..."
    ),
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
    h4("3. Return Periods"),
    textInput(
      "target_rperiods",
      "Target return periods (comma-separated):",
      value = "2, 5, 10, 20, 25, 50, 100, 150, 200"
    ),
    hr()
  )
}

# Create the data preview tab panel with statistics and table.
# Returns:
#   A Shiny tabPanel object displaying uploaded data summary,
#   interactive data table, and sample statistics.
create_data_preview_tab <- function() {
  tabPanel("Data Preview",
           h4("Sample Statistics"),
           uiOutput("sample_stats"),
           h4("Uploaded Data"),
           verbatimTextOutput("data_summary"),
           DTOutput("data_table"),
           hr()
  )
}

# Create the fitting results tab with parameter controls and plots.
# Returns:
#   A Shiny tabPanel object containing fitted parameters display,
#   manual parameter adjustment controls, confidence interval
#   computation options, main probability plot, GOF tests table,
#   quantiles table, and Excel report download button.
create_fitting_results_tab <- function() {
  tabPanel("Fitting Tool",
    fluidRow(
      column(12,
        div(
          style = "text-align: center; margin-bottom: 10px;",
          actionButton(
            "run_analysis", "Run Analysis",
            class = "btn-primary btn-md",
            icon = icon("play")
          )
        ),
        hr()
      )
    ),
    fluidRow(
      column(6,
        h4("Fitted Parameters (Initial Fit)"),
        verbatimTextOutput("initial_params"),
        hr(),
        h4("Manual Parameter Adjustment"),
        helpText(
          "Adjust parameters manually and click",
          "'Recompute' to update results"
        ),
        uiOutput("param_controls"),
        actionButton(
          "recompute",
          "Recompute with Manual Parameters",
          class = "btn-warning",
          icon = icon("calculator")
        ),
        hr(),
        h4("Confidence Intervals"),
        helpText(
          "Compute bootstrap confidence intervals for",
          "uncertainty quantification"
        ),
        fluidRow(
          column(
            4,
            numericInput(
              "ci_level", "Confidence Level:",
              value = 0.95, min = 0.5,
              max = 0.99, step = 0.01
            )
          ),
          column(
            4,
            numericInput(
              "n_bootstrap", "Bootstrap Iterations:",
              value = 100, min = 100,
              max = 10000, step = 100
            )
          ),
          column(
            4,
            checkboxInput(
              "parallel_bootstrap",
              "Use Parallel",
              value = FALSE
            )
          )
        ),
        actionButton(
          "compute_ci",
          "Compute Confidence Intervals",
          class = "btn-info",
          icon = icon("chart-line")
        ),
        hr()
      ),
      column(6,
        h4("Probability Plot (Value vs Return Period)"),
        downloadButton(
          "download_prob_rperiod",
          "Download Plot",
          class = "btn-sm"
        ),
        plotOutput("prob_plot_rperiod", height = "500px"),
        hr()
      )
    ),
    fluidRow(
      column(6,
        h4("Goodness-of-Fit Tests"),
        tableOutput("gof_tests")
      ),
      column(6,
        h4("Return Period Quantiles"),
        uiOutput("quantiles_table")
      )
    ),
    hr(),
    downloadButton(
      "download_report",
      "Download Excel Report",
      class = "btn-success"
    )
  )
}

# Create the diagnostic plots tab with multiple plot types.
# Returns:
#   A Shiny tabPanel object containing four diagnostic plots:
#   exceedance probability plot, Q-Q plot, histogram with fitted
#   PDF, and empirical vs fitted CDF. Each plot includes a
#   download button.
create_plots_tab <- function() {
  tabPanel("Auxiliary Plots",
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
    )
  )
}

# Create the main application UI layout.
# Returns:
#   A Shiny fluidPage object with sidebar layout containing all
#   tabs and components for the EVA Shiny application.
create_ui <- function() {
  fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("Extreme Value Analysis - Distribution Fitting"),
    sidebarLayout(
      create_sidebar_ui(),
      mainPanel(
        width = 9,
        tabsetPanel(
          id = "main_tabs",
          create_data_preview_tab(),
          create_fitting_results_tab(),
          create_plots_tab()
        )
      )
    )
  )
}
