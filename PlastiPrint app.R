# =============================================================================
# PlastiPrint — Shiny GUI for Microplastic Chemical Fingerprinting Workflow
# =============================================================================
# SPDX-License-Identifier: AGPL-3.0-only
#
# This Shiny application wraps the complete computational workflow from:
#   "Chemical Fingerprints of New vs. Weathered Microplastics:
#    A Machine Learning Approach"
#
# Each workflow step is exposed as a sequential tab with configurable parameters.
# =============================================================================

# ── 0.  Bootstrap: install required packages silently ──────────────────────────
local({
  cran_pkgs <- c(
    # Shiny + UI
    "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "DT",
    # Data manipulation
    "tidyverse", "readxl", "writexl", "data.table",
    # Visualization
    "crayon", "grid", "gridExtra", "ggplot2", "pheatmap",
    "viridis", "ggpubr", "scales",
    # Statistics
    "stats", "vegan", "cluster", "factoextra", "FactoMineR",
    "ggforce", "MASS", "moments",
    # Machine Learning
    "caret", "randomForest", "ranger", "e1071", "xgboost",
    # Missing data handling
    "missForest", "pROC", "VIM", "softImpute", "RANN",
    # Parallel processing
    "foreach", "doParallel", "parallel"
  )
  new <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
})

# ── 1.  Load libraries ────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(shinyWidgets)
  library(DT)
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(data.table)
  library(grid)
  library(gridExtra)
  library(stats)
  library(vegan)
  library(cluster)
  library(factoextra)
  library(FactoMineR)
  library(ggforce)
  library(MASS)
  library(caret)
  library(randomForest)
  library(ranger)
  library(e1071)
  library(xgboost)
  library(pheatmap)
  library(viridis)
  library(missForest)
  library(pROC)
  library(VIM)
  library(softImpute)
  library(ggpubr)
  library(RANN)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(ggplot2)
  library(moments)
  library(scales)
})

# ── 2.  Source helper functions ───────────────────────────────────────────────
helper_file <- file.path(getwd(), "Helper_Functions_for_Computational_Workflow_06Feb2026.R")
if (!file.exists(helper_file)) {
  # Also try one level up or in the app directory
  alt <- list.files(path = c(".", ".."), pattern = "Helper.*Functions.*\\.R$",
                    full.names = TRUE, recursive = FALSE)
  if (length(alt)) helper_file <- alt[1]
}
if (file.exists(helper_file)) {
  source(helper_file, local = FALSE)
} else {
  warning("Helper functions file not found. Place it in the app directory.")
}

# ── 2b.  Inline helper that is defined in the Rmd (not in helper file) ───────
print_filtering_stats <- function(results,
                                  meta_cols = c("Plastic_type", "technique",
                                                "Source", "Polymer")) {
  df_final <- results$filtered_data
  df_orig  <- results$original_data
  feats_final <- setdiff(names(df_final), meta_cols)
  feats_orig  <- setdiff(names(df_orig),  meta_cols)
  pct_features <- (length(feats_final) / length(feats_orig)) * 100
  lines <- character()
  lines <- c(lines, sprintf("Features retained: %d / %d (%.2f%%)",
                             length(feats_final), length(feats_orig), pct_features))
  if (length(feats_final) > 0) {
    vals_final <- df_final[, feats_final, drop = FALSE]
    pct_missing <- (sum(is.na(vals_final)) / prod(dim(vals_final))) * 100
    lines <- c(lines, sprintf("Missing values after reduction: %.2f%%", pct_missing))
  } else {
    lines <- c(lines, "Missing values after reduction: 0% (no features)")
    vals_final <- data.frame(0)
  }
  signal_final <- sum(vals_final, na.rm = TRUE)
  signal_orig  <- sum(df_orig[, feats_orig, drop = FALSE], na.rm = TRUE)
  pct_signal   <- if (signal_orig > 0) (signal_final / signal_orig) * 100 else 0
  lines <- c(lines, sprintf("Signal retained: %.2f%%", pct_signal))
  paste(lines, collapse = "\n")
}

# =============================================================================
#                               UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "PlastiPrint v1.0", titleWidth = 260),

  # ── Sidebar ──────────────────────────────────────────────────────────────────
 dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_tabs",
      menuItem("Welcome",           tabName = "welcome",  icon = icon("home")),
      menuItem("Step 1 — Import",   tabName = "step1",    icon = icon("file-import")),
      menuItem("Step 2 — Blank Sub",tabName = "step2",    icon = icon("filter")),
      menuItem("Step 3 — Source",   tabName = "step3",    icon = icon("tags")),
      menuItem("Step 4 — Metals",   tabName = "step4",    icon = icon("flask")),
      menuItem("Step 5 — Labels",   tabName = "step5",    icon = icon("link")),
      menuItem("Step 6 — Features", tabName = "step6",    icon = icon("cut")),
      menuItem("Step 7 — ML",       tabName = "step7",    icon = icon("brain")),
      menuItem("Step 8 — HCA",      tabName = "step8",    icon = icon("project-diagram"))
    ),
    hr(),
    div(style = "padding: 10px; font-size: 11px; color: #aaa;",
        "PlastiPrint wraps the complete",
        br(), "microplastic fingerprinting",
        br(), "computational workflow.")
  ),

  # ── Body ─────────────────────────────────────────────────────────────────────
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f7f9fc; }
      .box { border-top: 3px solid #3c8dbc; }
      .log-box { max-height: 350px; overflow-y: auto;
                 background: #1e1e1e; color: #d4d4d4;
                 padding: 12px; border-radius: 4px;
                 font-family: 'Fira Code', 'Consolas', monospace;
                 font-size: 12px; white-space: pre-wrap; }
      .nav-step-btn { margin-top: 15px; }
      .param-section { margin-bottom: 10px; }
      .small-label .control-label { font-size: 12px; }
    "))),

    tabItems(
      # ────────────────────── WELCOME ──────────────────────
      tabItem("welcome",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Welcome to PlastiPrint",
              h4("Chemical Fingerprints of New vs. Weathered Microplastics"),
              p("This GUI walks you through the complete computational",
                "workflow for microplastic chemical fingerprinting.",
                "Each tab on the left corresponds to a major step."),
              tags$ol(
                tags$li(tags$b("Step 1:"), " Import raw data (ATD-GC-MS / HPLC-TOF-MS) and group compounds"),
                tags$li(tags$b("Step 2:"), " Blank subtraction & negative-peak removal"),
                tags$li(tags$b("Step 3:"), " Assign source (Store-Bought / Environmental)"),
                tags$li(tags$b("Step 4:"), " Import ICP-MS trace-metal data"),
                tags$li(tags$b("Step 5:"), " Merge label information & create dataset combinations"),
                tags$li(tags$b("Step 6:"), " Feature reduction (Regular 80%, Modified 80%, In-House) & data fusion"),
                tags$li(tags$b("Step 7:"), " Hybrid pipeline ML optimization with ensemble feature selection"),
                tags$li(tags$b("Step 8:"), " Hierarchical Cluster Analysis")
              ),
              hr(),
              h4("Getting Started"),
              p("1. Set your ", tags$b("project root directory"), " below.",
                " This folder should contain the ", tags$code("Raw data/"),
                " subfolder and the helper-functions R file."),
              fluidRow(
                column(8, textInput("project_dir", "Project Root Directory:",
                                    value = getwd(), width = "100%")),
                column(4, actionButton("set_dir", "Set Directory",
                                       class = "btn-primary",
                                       style = "margin-top:25px;"))
              ),
              verbatimTextOutput("dir_status")
          )
        )
      ),

      # ────────────────────── STEP 1 ──────────────────────
      tabItem("step1",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 1: Data Import & Compound Grouping",
              h5("1.1 — Import Raw Data"),
              p("Reads CSV files from ", tags$code("Raw data/ATDGCMS"),
                " and XLS files from ", tags$code("Raw data/HPLCTOFMS"), "."),
              actionButton("run_step1_import", "Import Data",
                           class = "btn-primary", icon = icon("upload")),
              hr(),
              h5("1.2 — Compound Grouping Parameters"),
              div(class = "param-section",
                tags$b("ATD-GC-MS"),
                numericInput("gc_rtthres1", "RT Threshold:", 0.005,
                             min = 0.0001, step = 0.001, width = "60%"),
                numericInput("gc_mzthres", "m/z Threshold:", 0.05,
                             min = 0.0001, step = 0.01, width = "60%")
              ),
              div(class = "param-section",
                tags$b("HPLC-TOF-MS"),
                numericInput("hplc_split_time", "Split Time (min):", 6,
                             min = 1, step = 0.5, width = "60%"),
                numericInput("hplc_rtthres1", "RT Threshold 1 (≤split):", 0.4,
                             min = 0.01, step = 0.05, width = "60%"),
                numericInput("hplc_rtthres2", "RT Threshold 2 (>split):", 0.28,
                             min = 0.01, step = 0.05, width = "60%"),
                numericInput("hplc_mzthres", "m/z Threshold:", 0.0005,
                             min = 0.00001, step = 0.0001, width = "60%")
              ),
              actionButton("run_step1_group", "Group Compounds",
                           class = "btn-success", icon = icon("object-group")),
              hr(),
              h5("1.3 — Benchmark Removal"),
              p("Removes known benchmark/internal-standard compounds."),
              actionButton("run_step1_bench", "Remove Benchmarks",
                           class = "btn-warning", icon = icon("trash")),
              hr(),
              actionButton("go_step2", "Next → Step 2",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output / Preview",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step1_log"))),
                tabPanel("GC Preview",  DTOutput("step1_gc_table")),
                tabPanel("HPLC Preview", DTOutput("step1_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 2 ──────────────────────
      tabItem("step2",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 2: Blank Subtraction",
              p("For each feature, the blank mean is subtracted from sample",
                "areas/heights. Values below 3× blank SD are removed."),
              numericInput("blank_sd_mult", "Blank SD Multiplier:", 3,
                           min = 1, max = 10, step = 0.5),
              actionButton("run_step2", "Run Blank Subtraction",
                           class = "btn-primary", icon = icon("filter")),
              hr(),
              actionButton("go_step3", "Next → Step 3",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step2_log"))),
                tabPanel("GC After Blank Sub", DTOutput("step2_gc_table")),
                tabPanel("HPLC After Blank Sub", DTOutput("step2_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 3 ──────────────────────
      tabItem("step3",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 3: Assign Source & Rename",
              p("Labels each sample as ", tags$b("Store-Bought"),
                " or ", tags$b("Environmental"),
                " and assigns an instrument technique tag."),
              p("Environmental samples contain 'USE' in their filename."),
              actionButton("run_step3", "Run Step 3",
                           class = "btn-primary", icon = icon("tags")),
              hr(),
              actionButton("go_step4", "Next → Step 4",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step3_log"))),
                tabPanel("GC Data", DTOutput("step3_gc_table")),
                tabPanel("HPLC Data", DTOutput("step3_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 4 ──────────────────────
      tabItem("step4",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 4: Import Trace-Metal (ICP-MS) Data",
              p("Reads ICP-MS data from Excel files, selects best",
                "analytical mode per element, removes seawater metals."),
              textInput("metals_remove", "Seawater metals to discard:",
                        "Na, Mg, K, Ca"),
              actionButton("run_step4", "Import Metals",
                           class = "btn-primary", icon = icon("flask")),
              hr(),
              actionButton("go_step5", "Next → Step 5",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step4_log"))),
                tabPanel("ICP-MS Data", DTOutput("step4_icp_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 5 ──────────────────────
      tabItem("step5",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 5: Combine Label Info & Create Datasets",
              p("Merges sample label info (Plastic_type, Subcategory,",
                "Polymer) and creates all single- and multi-instrument",
                "dataset combinations."),
              fileInput("label_file", "Sample Label Excel File (.xlsx):",
                        accept = c(".xlsx", ".xls")),
              p(tags$small("Leave empty to use default file in Raw data/.")),
              actionButton("run_step5", "Run Step 5",
                           class = "btn-primary", icon = icon("link")),
              hr(),
              actionButton("go_step6", "Next → Step 6",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step5_log"))),
                tabPanel("Datasets Summary", DTOutput("step5_summary_table")),
                tabPanel("GC Labeled", DTOutput("step5_gc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 6 ──────────────────────
      tabItem("step6",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 6: Feature Reduction & Data Fusion",
              h5("6.1–6.3  Feature Reduction"),
              selectInput("feat_dataset", "Dataset to filter:",
                          choices = c("gc", "hplc", "icp",
                                      "gc_hplc", "gc_icp",
                                      "hplc_icp", "gc_hplc_icp"),
                          selected = "gc"),
              radioButtons("feat_mode", "Filtering method:",
                           choices = c("Regular 80% rule" = "regular80",
                                       "Modified 80% rule" = "modified80",
                                       "In-house comprehensive" = "inhouse"),
                           selected = "inhouse"),
              actionButton("run_step6_filter", "Run Feature Reduction",
                           class = "btn-primary", icon = icon("cut")),
              hr(),
              h5("6.4  Data Fusion"),
              p("Applies in-house filtering to GC, HPLC, ICP",
                "individually, then fuses multi-instrument datasets."),
              actionButton("run_step6_fusion", "Run Data Fusion",
                           class = "btn-success", icon = icon("compress")),
              hr(),
              actionButton("go_step7", "Next → Step 7",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step6_log"))),
                tabPanel("Filtering Stats",
                         verbatimTextOutput("step6_stats")),
                tabPanel("Filtered Data Preview", DTOutput("step6_table")),
                tabPanel("Fusion Summary", DTOutput("step6_fusion_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 7 ──────────────────────
      tabItem("step7",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 7: Hybrid Pipeline ML Optimization",
              h5("Train / Test Strategy"),
              radioButtons("step7_option", "Training approach:",
                choices = c(
                  "Store-Bought → Predict Environmental" = "sb_only",
                  "Both (SB+ENV) → Predict Environmental" = "sb_env"
                ), selected = "sb_only"),
              h5("Datasets to analyze"),
              checkboxGroupInput("step7_datasets", NULL,
                choices = c("gc", "gc_icp"),
                selected = c("gc", "gc_icp")),
              hr(),
              h5("Pipeline Search Space"),
              checkboxGroupInput("step7_impute", "Imputation methods:",
                choices  = c("half_min", "median", "knn"),
                selected = c("half_min", "median", "knn"), inline = TRUE),
              checkboxGroupInput("step7_norm", "Normalization methods:",
                choices  = c("none", "log", "log10", "zscore", "pareto"),
                selected = c("none", "log", "log10", "zscore", "pareto"),
                inline = TRUE),
              checkboxGroupInput("step7_algo", "ML algorithms:",
                choices  = c("ranger", "svmRadial", "xgbTree"),
                selected = c("ranger", "svmRadial", "xgbTree"), inline = TRUE),
              hr(),
              h5("Feature Selection Thresholds"),
              sliderInput("step7_top_pct", "Top % features per pipeline:",
                          min = 0.05, max = 0.50, value = 0.20, step = 0.05),
              sliderInput("step7_stab", "Stability threshold:",
                          min = 0.30, max = 0.90, value = 0.60, step = 0.05),
              sliderInput("step7_rank", "Mean rank percentile:",
                          min = 0.10, max = 0.50, value = 0.30, step = 0.05),
              hr(),
              h5("Other Settings"),
              numericInput("step7_cv", "CV folds:", 5, min = 2, max = 10),
              numericInput("step7_topk", "Top-k pipelines:", 5, min = 2, max = 20),
              numericInput("step7_nperm", "Permutation reps:", 10, min = 5, max = 50),
              numericInput("step7_seed", "Random seed:", 123, min = 1),
              checkboxInput("step7_parallel", "Use parallel processing", FALSE),
              hr(),
              actionButton("run_step7", "Run Hybrid Analysis",
                           class = "btn-danger btn-lg",
                           icon = icon("play"),
                           style = "width:100%;"),
              p(tags$small("⚠ This step may take several minutes.")),
              hr(),
              actionButton("go_step8", "Next → Step 8",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Results",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step7_log"))),
                tabPanel("Best Pipeline", verbatimTextOutput("step7_best")),
                tabPanel("Metrics", DTOutput("step7_metrics_table")),
                tabPanel("Confusion Matrix", plotOutput("step7_cm_plot",
                                                        height = "500px")),
                tabPanel("Feature Tiers", plotOutput("step7_tier_plot",
                                                     height = "600px")),
                tabPanel("Top Features", DTOutput("step7_features_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 8 ──────────────────────
      tabItem("step8",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 8: Hierarchical Cluster Analysis",
              selectInput("hca_source", "Use results from:",
                choices = c("Step 7 Option 1 (SB train)" = "sb",
                            "Step 7 Option 2 (SB+ENV train)" = "sb_env"),
                selected = "sb"),
              selectInput("hca_dataset", "Dataset:", choices = c("gc", "gc_icp"),
                          selected = "gc"),
              selectInput("hca_dist", "Distance method:",
                choices = c("manhattan", "euclidean", "bray"),
                selected = "manhattan"),
              selectInput("hca_link", "Linkage method:",
                choices = c("average", "complete", "ward.D2", "single"),
                selected = "average"),
              actionButton("run_step8", "Run HCA",
                           class = "btn-primary", icon = icon("project-diagram")),
              hr(),
              h5("Export"),
              downloadButton("download_hca", "Download Dendrogram (PDF)")
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Dendrogram",
              tabsetPanel(
                tabPanel("Plot",
                         plotOutput("step8_dendro", height = "600px")),
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step8_log"))),
                tabPanel("CCC",
                         verbatimTextOutput("step8_ccc"))
              )
          )
        )
      )
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# =============================================================================
#                              SERVER
# =============================================================================
server <- function(input, output, session) {

  # ── Reactive storage for intermediate data ──────────────────────────────────
  rv <- reactiveValues(
    # Step 1
    atdgcms_raw     = NULL,
    atdgcms_blank   = NULL,
    atdgcms_step1   = NULL,
    all_hplc        = NULL,
    atdgcms_grouped = NULL,
    hplc_grouped    = NULL,
    # Step 2
    atdgcms_step2   = NULL,
    hplc_step2      = NULL,
    # Step 3
    gc              = NULL,
    hplc            = NULL,
    # Step 4
    icp             = NULL,
    # Step 5
    sampinfo        = NULL,
    gc_labeled      = NULL,
    hplc_labeled    = NULL,
    icp_labeled     = NULL,
    gc_hplc         = NULL,
    gc_icp          = NULL,
    hplc_icp        = NULL,
    gc_hplc_icp     = NULL,
    # Step 6
    filter_results  = NULL,
    gc_clean        = NULL,
    hplc_clean      = NULL,
    icp_clean       = NULL,
    gc_icp_clean    = NULL,
    gc_hplc_clean   = NULL,
    hplc_icp_clean  = NULL,
    gc_hplc_icp_clean = NULL,
    # Step 7
    hybrid_results_sb     = NULL,
    hybrid_results_sb_env = NULL,
    # Step 8
    hca_result      = NULL,

    # Logs
    log_step1 = "Awaiting input...",
    log_step2 = "Run Step 1 first.",
    log_step3 = "Run Step 2 first.",
    log_step4 = "Awaiting input...",
    log_step5 = "Run Steps 1-4 first.",
    log_step6 = "Run Step 5 first.",
    log_step7 = "Run Step 6 first.",
    log_step8 = "Run Step 7 first.",
    step6_stats_text = ""
  )

  # ── Helper: append to a log ─────────────────────────────────────────────────
  append_log <- function(key, msg) {
    old <- isolate(rv[[key]])
    if (is.null(old) || old == "Awaiting input..." ||
        grepl("^Run Step", old)) {
      rv[[key]] <- msg
    } else {
      rv[[key]] <- paste0(old, "\n", msg)
    }
  }

  # ── Navigation helpers ──────────────────────────────────────────────────────
  observeEvent(input$go_step2, updateTabItems(session, "main_tabs", "step2"))
  observeEvent(input$go_step3, updateTabItems(session, "main_tabs", "step3"))
  observeEvent(input$go_step4, updateTabItems(session, "main_tabs", "step4"))
  observeEvent(input$go_step5, updateTabItems(session, "main_tabs", "step5"))
  observeEvent(input$go_step6, updateTabItems(session, "main_tabs", "step6"))
  observeEvent(input$go_step7, updateTabItems(session, "main_tabs", "step7"))
  observeEvent(input$go_step8, updateTabItems(session, "main_tabs", "step8"))

  # ── Set directory ───────────────────────────────────────────────────────────
  observeEvent(input$set_dir, {
    d <- trimws(input$project_dir)
    if (dir.exists(d)) {
      setwd(d)
      output$dir_status <- renderText(paste("✓ Working directory set to:", getwd()))
    } else {
      output$dir_status <- renderText(paste("✗ Directory not found:", d))
    }
  })

  # ===========================================================================
  #  STEP 1  —  Import, Group, Benchmark
  # ===========================================================================

  # ── 1.1: Import raw data ───────────────────────────────────────────────────
  observeEvent(input$run_step1_import, {
    rv$log_step1 <- ""
    tryCatch({
      # --- ATD-GC-MS ---
      append_log("log_step1", "── Importing ATD-GC-MS data ──")
      path_atdgcms <- "./Raw data/ATDGCMS"
      if (!dir.exists(path_atdgcms)) stop(paste("Directory not found:", path_atdgcms))

      atdgcms_list <- list.files(path_atdgcms, pattern = "*.csv", full.names = TRUE) %>%
        .[!str_detect(., "2022")] %>%
        .[!str_detect(., "USE_(01|02|05|06|09|10|11|12|13|15|16|17|18|19|20)")] %>%
        .[!str_detect(., "USSB_(01|08)")]
      blank_list <- list.files(path_atdgcms, pattern = "*.csv", full.names = TRUE) %>%
        .[str_detect(., "2022-")]

      read_and_add_filename <- function(file) {
        data.table::fread(file, showProgress = FALSE, data.table = FALSE,
                          check.names = TRUE) %>%
          mutate(File = basename(file))
      }

      atdgcms_raw <- data.table::rbindlist(
        lapply(atdgcms_list, read_and_add_filename),
        use.names = TRUE, fill = TRUE) %>%
        as.data.frame() %>%
        dplyr::select(-any_of(c("Start","End","Width","Base.Peak",
                                "Cpd","Label","Height","Ions"))) %>%
        dplyr::mutate(File = gsub("_", "-", File)) %>%
        dplyr::mutate(File = gsub(".csv", "", File)) %>%
        mutate(type = "Sample")

      atdgcms_blank <- data.table::rbindlist(
        lapply(blank_list, function(path) {
          data.table::fread(path, showProgress = FALSE, check.names = TRUE)
        }), use.names = TRUE, fill = TRUE) %>%
        as.data.frame() %>%
        dplyr::select(any_of(c("Area","File","m.z","RT"))) %>%
        mutate(type = "Blanks")

      atdgcms_step1 <- bind_rows(atdgcms_raw, atdgcms_blank) %>%
        arrange(RT) %>%
        mutate(simplified_file = map_chr(File, ~ paste(
          str_split_1(.x, pattern = "-")[c(4, 5)], collapse = "-")))

      rv$atdgcms_raw   <- atdgcms_raw
      rv$atdgcms_blank <- atdgcms_blank
      rv$atdgcms_step1 <- atdgcms_step1
      append_log("log_step1", sprintf("  GC samples: %d rows from %d files",
                                      nrow(atdgcms_raw), length(atdgcms_list)))
      append_log("log_step1", sprintf("  GC blanks : %d rows from %d files",
                                      nrow(atdgcms_blank), length(blank_list)))

      # --- HPLC-TOF-MS ---
      append_log("log_step1", "\n── Importing HPLC-TOF-MS data ──")
      path_batch0 <- "./Raw data/HPLCTOFMS/EF_Non-target data_Batch 0"
      path_batch1 <- "./Raw data/HPLCTOFMS/EF_Non-target data_Batch 1"
      path_batch2 <- "./Raw data/HPLCTOFMS/EF_Non-target data_Batch 2"

      read_hplc_batch_gui <- function(dir_path, batch_num) {
        if (!dir.exists(dir_path)) {
          warning(paste("Directory not found:", dir_path))
          return(data.frame())
        }
        files_all <- list.files(dir_path, pattern = "*.xls", full.names = TRUE)
        if (length(files_all) == 0) return(data.frame())

        read_xls_fast <- function(path) readxl::read_xls(path, skip = 1) %>% as.data.frame()

        combined <- data.table::rbindlist(
          lapply(files_all, read_xls_fast),
          use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          dplyr::select(m.z = `m/z`, RT, Height, File) %>%
          mutate(File = gsub("_", "-", File),
                 type = ifelse(str_detect(File, "Blank"), "Blanks", "Sample"),
                 Day = case_when(str_detect(File, "Day0") ~ "0",
                                 str_detect(File, "Day2") ~ "2",
                                 str_detect(File, "Day15") ~ "15",
                                 TRUE ~ "34"),
                 Replicate = ifelse(str_detect(File, "R1"), "R1", "R2"),
                 Batch_number = as.character(batch_num))
        return(combined)
      }

      hplc_b0 <- read_hplc_batch_gui(path_batch0, 0) %>%
        mutate(simplified_file = map_chr(File, ~ paste(str_split_1(.x, "-")[1:2], collapse = "-")))
      hplc_b1 <- read_hplc_batch_gui(path_batch1, 1) %>%
        mutate(simplified_file = map_chr(File, ~ paste(str_split_1(.x, "-")[2:3], collapse = "-")))
      hplc_b2 <- read_hplc_batch_gui(path_batch2, 2) %>%
        mutate(simplified_file = map_chr(File, ~ paste(str_split_1(.x, "-")[2:3], collapse = "-")))

      all_hplc <- data.table::rbindlist(list(hplc_b0, hplc_b1, hplc_b2),
                                        use.names = TRUE, fill = TRUE) %>%
        as.data.frame() %>% arrange(RT)
      rv$all_hplc <- all_hplc
      append_log("log_step1", sprintf("  HPLC total: %d rows", nrow(all_hplc)))
      append_log("log_step1", "\n✓ Data import complete.")
    },
    error = function(e) {
      append_log("log_step1", paste("✗ ERROR:", conditionMessage(e)))
    })
  })

  # ── 1.2: Compound grouping ─────────────────────────────────────────────────
  observeEvent(input$run_step1_group, {
    req(rv$atdgcms_step1, rv$all_hplc)
    tryCatch({
      append_log("log_step1", "\n── Grouping compounds ──")
      rv$atdgcms_grouped <- grouping_comp(
        rv$atdgcms_step1,
        rtthres1 = input$gc_rtthres1,
        mzthres  = input$gc_mzthres,
        type     = "ATDGCMS")
      append_log("log_step1", sprintf("  GC grouped: %d unique features",
                                      n_distinct(rv$atdgcms_grouped$Feature, na.rm = TRUE)))

      rv$hplc_grouped <- grouping_comp(
        rv$all_hplc,
        split_time = input$hplc_split_time,
        rtthres1   = input$hplc_rtthres1,
        rtthres2   = input$hplc_rtthres2,
        mzthres    = input$hplc_mzthres,
        type       = "HPLCTOFMS")
      append_log("log_step1", sprintf("  HPLC grouped: %d unique features",
                                      n_distinct(rv$hplc_grouped$Feature, na.rm = TRUE)))
      append_log("log_step1", "✓ Grouping complete.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  # ── 1.3: Benchmark removal ─────────────────────────────────────────────────
  observeEvent(input$run_step1_bench, {
    req(rv$atdgcms_grouped, rv$hplc_grouped)
    tryCatch({
      append_log("log_step1", "\n── Removing benchmarks ──")

      gc_g <- rv$atdgcms_grouped
      n_before_gc <- nrow(gc_g)
      idx <- which(gc_g$m.z >= 98 & gc_g$m.z <= 98.5 &
                   gc_g$RT > 13.4 & gc_g$RT < 13.5)
      if (length(idx) > 0) gc_g <- gc_g[-idx, ]
      rv$atdgcms_grouped <- gc_g
      append_log("log_step1", sprintf("  GC: removed %d benchmark rows", n_before_gc - nrow(gc_g)))

      hplc_g <- rv$hplc_grouped
      n_before_hplc <- nrow(hplc_g)
      hplc_benchmark_mz <- c(156.0957, 198.1723, 246.2568, 305.1956, 342.1722)
      hplc_idx <- c()
      for (mz in hplc_benchmark_mz) {
        idx <- which(abs(hplc_g$m.z - mz) <= 0.0005)
        hplc_idx <- c(hplc_idx, idx)
      }
      if (length(hplc_idx) > 0) hplc_g <- hplc_g[-hplc_idx, ]
      rv$hplc_grouped <- hplc_g
      append_log("log_step1", sprintf("  HPLC: removed %d benchmark rows",
                                      n_before_hplc - nrow(hplc_g)))
      append_log("log_step1", "✓ Benchmarks removed.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  # Step 1 outputs
  output$step1_log <- renderText(rv$log_step1)
  output$step1_gc_table <- renderDT({
    req(rv$atdgcms_grouped)
    datatable(head(rv$atdgcms_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step1_hplc_table <- renderDT({
    req(rv$hplc_grouped)
    datatable(head(rv$hplc_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 2  —  Blank Subtraction
  # ===========================================================================
  observeEvent(input$run_step2, {
    req(rv$atdgcms_grouped, rv$hplc_grouped)
    rv$log_step2 <- ""
    tryCatch({
      sd_mult <- input$blank_sd_mult
      append_log("log_step2", sprintf("── Blank subtraction (SD multiplier = %.1f) ──", sd_mult))

      # GC
      gc_adj <- rv$atdgcms_grouped %>%
        dplyr::group_by(Feature) %>%
        dplyr::mutate(
          has_blank = any(type == "Blanks"),
          blank_mean = mean(Area[type == "Blanks"], na.rm = TRUE),
          blank_sd   = sd(Area[type == "Blanks"], na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::filter(type != "Blanks") %>%
        dplyr::mutate(
          Area = if_else(has_blank, Area - blank_mean, Area),
          threshold = if_else(has_blank, sd_mult * blank_sd, 0)) %>%
        dplyr::filter(Area > 0, Area > threshold) %>%
        dplyr::select(-has_blank, -blank_mean, -blank_sd, -threshold)
      rv$atdgcms_step2 <- gc_adj
      append_log("log_step2", sprintf("  GC after blank sub: %d rows", nrow(gc_adj)))

      # HPLC
      hplc_adj <- rv$hplc_grouped %>%
        dplyr::group_by(Batch_number, Feature) %>%
        dplyr::mutate(
          has_blank = any(type == "Blanks"),
          blank_mean = mean(Height[type == "Blanks"], na.rm = TRUE),
          blank_sd   = sd(Height[type == "Blanks"], na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::filter(type != "Blanks") %>%
        dplyr::mutate(
          Height = if_else(has_blank, Height - blank_mean, Height),
          threshold = if_else(has_blank, sd_mult * blank_sd, 0)) %>%
        dplyr::filter(Height > 0, Height > threshold) %>%
        dplyr::select(-has_blank, -blank_mean, -blank_sd, -threshold)
      rv$hplc_step2 <- hplc_adj
      append_log("log_step2", sprintf("  HPLC after blank sub: %d rows", nrow(hplc_adj)))
      append_log("log_step2", "✓ Blank subtraction complete.")
    },
    error = function(e) append_log("log_step2", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step2_log <- renderText(rv$log_step2)
  output$step2_gc_table <- renderDT({
    req(rv$atdgcms_step2)
    datatable(head(rv$atdgcms_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step2_hplc_table <- renderDT({
    req(rv$hplc_step2)
    datatable(head(rv$hplc_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 3  —  Source assignment
  # ===========================================================================
  observeEvent(input$run_step3, {
    req(rv$atdgcms_step2, rv$hplc_step2)
    rv$log_step3 <- ""
    tryCatch({
      append_log("log_step3", "── Assigning Source & Technique ──")

      gc <- rv$atdgcms_step2 %>%
        dplyr::filter(!is.na(Feature)) %>%
        dplyr::select(File, Feature, m.z, RT, Area) %>%
        mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
        dplyr::rename(Values = Area) %>%
        mutate(technique = "GC")
      rv$gc <- gc

      hplc <- rv$hplc_step2 %>%
        dplyr::filter(!is.na(Feature)) %>%
        dplyr::select(File, Feature, m.z, RT, Height) %>%
        dplyr::mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
        dplyr::rename(Values = Height) %>%
        mutate(technique = "HPLC")
      rv$hplc <- hplc

      append_log("log_step3", sprintf("  GC:   %d rows  |  SB: %d  ENV: %d",
        nrow(gc), sum(gc$Source == "Store-Bought"), sum(gc$Source == "Environmental")))
      append_log("log_step3", sprintf("  HPLC: %d rows  |  SB: %d  ENV: %d",
        nrow(hplc), sum(hplc$Source == "Store-Bought"), sum(hplc$Source == "Environmental")))
      append_log("log_step3", "✓ Step 3 complete.")
    },
    error = function(e) append_log("log_step3", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step3_log <- renderText(rv$log_step3)
  output$step3_gc_table <- renderDT({
    req(rv$gc)
    datatable(head(rv$gc, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step3_hplc_table <- renderDT({
    req(rv$hplc)
    datatable(head(rv$hplc, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 4  —  ICP-MS Import
  # ===========================================================================
  observeEvent(input$run_step4, {
    rv$log_step4 <- ""
    tryCatch({
      append_log("log_step4", "── Importing ICP-MS trace-metal data ──")

      read_trace_metal <- function(path) {
        if (!file.exists(path)) stop(paste("File missing:", path))
        readxl::read_excel(path) %>% column_to_rownames(var = "File")
      }
      remove_H2_HMI_modus <- function(icpms_df) {
        strings <- colnames(icpms_df)
        pattern <- "^(.*\\[ )([A-Za-z0-9 ]+)( \\])$"
        keep_cols <- strings[!str_detect(str_match(strings, pattern)[,3], "H2 HMI")]
        icpms_df %>% dplyr::select(all_of(keep_cols))
      }
      get_final_metal_names <- function(recovery_uncertainty, target_cols) {
        recovery_uncertainty <- recovery_uncertainty %>% filter(`Element+Modus` %in% target_cols)
        first_numbers <- sub(" .*", "", recovery_uncertainty$`Element+Modus`)
        dupes <- first_numbers[duplicated(first_numbers)] %>% unique()
        kept_names <- c()
        for (num in unique(first_numbers)) {
          rows <- recovery_uncertainty[first_numbers == num, ]
          if (nrow(rows) == 1) { kept_names <- c(kept_names, rows$`Element+Modus`)
          } else {
            best_row <- rows[which.max(rows$`Recovery [%]`), ]
            if (nrow(best_row) == 0) best_row <- rows[which.min(rows$`U(k=2)`), ]
            kept_names <- c(kept_names, best_row$`Element+Modus`)
          }
        }
        return(kept_names)
      }

      icpms_r1 <- read_trace_metal("./Raw data/ICPMS_Trace metal/icpms_round1_rawdata_removal USE-01-rep1 and USE-03 (only2observations).xlsx") %>% remove_H2_HMI_modus()
      icpms_r2b1 <- read_trace_metal("./Raw data/ICPMS_Trace metal/icpms_round2_batch1_rawdata.xlsx") %>% remove_H2_HMI_modus()
      icpms_r2b2 <- read_trace_metal("./Raw data/ICPMS_Trace metal/icpms_round2_batch2_rawdata.xlsx") %>% remove_H2_HMI_modus()

      rec_unc_1 <- readxl::read_excel("./Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round1.xlsx")
      rec_unc_2_1 <- readxl::read_excel("./Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round2_batch1.xlsx")
      rec_unc_2_2 <- readxl::read_excel("./Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round2_batch2.xlsx")

      fn1 <- get_final_metal_names(rec_unc_1, colnames(icpms_r1))
      fn2_1 <- get_final_metal_names(rec_unc_2_1, colnames(icpms_r2b1))
      fn2_2 <- get_final_metal_names(rec_unc_2_2, colnames(icpms_r2b2))

      icpms_r1   <- icpms_r1 %>% dplyr::select(all_of(fn1))
      icpms_r2b1 <- icpms_r2b1 %>% dplyr::select(all_of(fn2_1))
      icpms_r2b2 <- icpms_r2b2 %>% dplyr::select(all_of(fn2_2))

      common_cols <- Reduce(intersect, list(names(icpms_r1), names(icpms_r2b1), names(icpms_r2b2)))
      combined_icpms <- rbind(icpms_r1[common_cols], icpms_r2b1[common_cols], icpms_r2b2[common_cols])

      metals_remove <- trimws(unlist(strsplit(input$metals_remove, ",")))
      cols_to_remove <- unique(unlist(sapply(metals_remove, function(x) grep(x, colnames(combined_icpms)))))
      if (length(cols_to_remove) > 0) combined_icpms <- combined_icpms[, -cols_to_remove]

      icp <- combined_icpms %>%
        rownames_to_column(var = "File") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Feature", values_to = "Values") %>%
        mutate(technique = "ICP",
               Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought"),
               Values = as.numeric(Values))
      rv$icp <- icp

      append_log("log_step4", sprintf("  ICP-MS: %d rows, %d elements",
                                      nrow(icp), n_distinct(icp$Feature)))
      append_log("log_step4", "✓ ICP-MS import complete.")
    },
    error = function(e) append_log("log_step4", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step4_log <- renderText(rv$log_step4)
  output$step4_icp_table <- renderDT({
    req(rv$icp)
    datatable(head(rv$icp, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 5  —  Label merge & Combinations
  # ===========================================================================
  observeEvent(input$run_step5, {
    req(rv$gc, rv$hplc, rv$icp)
    rv$log_step5 <- ""
    tryCatch({
      append_log("log_step5", "── Merging label information ──")

      label_path <- if (!is.null(input$label_file)) {
        input$label_file$datapath
      } else {
        "./Raw data/Sample Labelling_all data_GC+HPLC+ICP_30Sept2025.xlsx"
      }
      if (!file.exists(label_path)) stop("Label file not found: ", label_path)

      sampinfo <- readxl::read_excel(label_path) %>%
        dplyr::select(c(File, Plastic_type, Subcategory, Polymer))
      rv$sampinfo <- sampinfo

      gc   <- left_join(rv$gc,   sampinfo, by = "File")
      hplc <- left_join(rv$hplc, sampinfo, by = "File")
      icp  <- left_join(rv$icp,  sampinfo, by = "File")

      rv$gc_labeled   <- gc
      rv$hplc_labeled <- hplc
      rv$icp_labeled  <- icp

      gc_hplc     <- bind_rows(gc, hplc)
      gc_icp      <- bind_rows(icp, gc %>% dplyr::select(-c("m.z", "RT")))
      hplc_icp    <- bind_rows(icp, hplc %>% dplyr::select(-c("m.z", "RT")))
      gc_hplc_icp <- bind_rows(icp, gc_hplc %>% dplyr::select(-c("m.z", "RT")))

      rv$gc_hplc     <- gc_hplc
      rv$gc_icp      <- gc_icp
      rv$hplc_icp    <- hplc_icp
      rv$gc_hplc_icp <- gc_hplc_icp

      append_log("log_step5", sprintf("  GC:   %d rows  (%d plastic types)",
                                      nrow(gc), n_distinct(gc$Plastic_type, na.rm = TRUE)))
      append_log("log_step5", sprintf("  HPLC: %d rows  (%d plastic types)",
                                      nrow(hplc), n_distinct(hplc$Plastic_type, na.rm = TRUE)))
      append_log("log_step5", sprintf("  ICP:  %d rows  (%d plastic types)",
                                      nrow(icp), n_distinct(icp$Plastic_type, na.rm = TRUE)))
      append_log("log_step5", "  Combinations created: gc_hplc, gc_icp, hplc_icp, gc_hplc_icp")
      append_log("log_step5", "✓ Step 5 complete.")
    },
    error = function(e) append_log("log_step5", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step5_log <- renderText(rv$log_step5)
  output$step5_summary_table <- renderDT({
    req(rv$gc_labeled)
    summary_df <- data.frame(
      Dataset = c("gc", "hplc", "icp", "gc_hplc", "gc_icp", "hplc_icp", "gc_hplc_icp"),
      Rows = c(nrow(rv$gc_labeled), nrow(rv$hplc_labeled), nrow(rv$icp_labeled),
               nrow(rv$gc_hplc), nrow(rv$gc_icp), nrow(rv$hplc_icp), nrow(rv$gc_hplc_icp)),
      Features = c(n_distinct(rv$gc_labeled$Feature), n_distinct(rv$hplc_labeled$Feature),
                   n_distinct(rv$icp_labeled$Feature), n_distinct(rv$gc_hplc$Feature),
                   n_distinct(rv$gc_icp$Feature), n_distinct(rv$hplc_icp$Feature),
                   n_distinct(rv$gc_hplc_icp$Feature))
    )
    datatable(summary_df, options = list(dom = 't'))
  })
  output$step5_gc_table <- renderDT({
    req(rv$gc_labeled)
    datatable(head(rv$gc_labeled, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 6  —  Feature Reduction & Data Fusion
  # ===========================================================================

  # Helper to get dataset by name
  get_dataset <- function(name) {
    switch(name,
      gc          = rv$gc_labeled,
      hplc        = rv$hplc_labeled,
      icp         = rv$icp_labeled,
      gc_hplc     = rv$gc_hplc,
      gc_icp      = rv$gc_icp,
      hplc_icp    = rv$hplc_icp,
      gc_hplc_icp = rv$gc_hplc_icp,
      NULL)
  }

  # ── 6.1–6.3: Feature reduction ─────────────────────────────────────────────
  observeEvent(input$run_step6_filter, {
    rv$log_step6 <- ""
    tryCatch({
      ds_name <- input$feat_dataset
      ds <- get_dataset(ds_name)
      req(ds)
      mode <- input$feat_mode

      append_log("log_step6", sprintf("── Feature reduction: %s [%s] ──", ds_name, mode))

      results <- shared_feature_filtering(ds, feat_reduc_mode = mode)
      rv$filter_results <- results

      stats_text <- print_filtering_stats(results)
      rv$step6_stats_text <- stats_text
      append_log("log_step6", stats_text)
      append_log("log_step6", "✓ Feature reduction complete.")
    },
    error = function(e) append_log("log_step6", paste("✗ ERROR:", conditionMessage(e))))
  })

  # ── 6.4: Data fusion ───────────────────────────────────────────────────────
  observeEvent(input$run_step6_fusion, {
    req(rv$gc_labeled, rv$hplc_labeled, rv$icp_labeled)
    tryCatch({
      append_log("log_step6", "\n── Running in-house filter + data fusion ──")

      gc_clean   <- shared_feature_filtering(rv$gc_labeled,   feat_reduc_mode = "inhouse")$filtered_data
      hplc_clean <- shared_feature_filtering(rv$hplc_labeled, feat_reduc_mode = "inhouse")$filtered_data
      icp_clean  <- shared_feature_filtering(rv$icp_labeled,  feat_reduc_mode = "inhouse")$filtered_data
      rv$gc_clean   <- gc_clean
      rv$hplc_clean <- hplc_clean
      rv$icp_clean  <- icp_clean

      rv$gc_hplc_clean     <- data_fusion(rv$gc_hplc,     source_feats = list(gc_clean, hplc_clean))
      rv$gc_icp_clean      <- data_fusion(rv$gc_icp,      source_feats = list(gc_clean, icp_clean))
      rv$hplc_icp_clean    <- data_fusion(rv$hplc_icp,    source_feats = list(hplc_clean, icp_clean))
      rv$gc_hplc_icp_clean <- data_fusion(rv$gc_hplc_icp, source_feats = list(gc_clean, hplc_clean, icp_clean))

      meta <- c("Plastic_type", "technique", "Source", "Polymer")
      fusion_summary <- data.frame(
        Dataset = c("gc", "hplc", "icp", "gc_hplc", "gc_icp", "hplc_icp", "gc_hplc_icp"),
        Samples = c(nrow(gc_clean), nrow(hplc_clean), nrow(icp_clean),
                    nrow(rv$gc_hplc_clean), nrow(rv$gc_icp_clean),
                    nrow(rv$hplc_icp_clean), nrow(rv$gc_hplc_icp_clean)),
        Features = c(
          length(setdiff(names(gc_clean), meta)),
          length(setdiff(names(hplc_clean), meta)),
          length(setdiff(names(icp_clean), meta)),
          length(setdiff(names(rv$gc_hplc_clean), meta)),
          length(setdiff(names(rv$gc_icp_clean), meta)),
          length(setdiff(names(rv$hplc_icp_clean), meta)),
          length(setdiff(names(rv$gc_hplc_icp_clean), meta))
        )
      )
      rv$fusion_summary <- fusion_summary
      append_log("log_step6", "✓ Data fusion complete. See 'Fusion Summary' tab.")
    },
    error = function(e) append_log("log_step6", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step6_log <- renderText(rv$log_step6)
  output$step6_stats <- renderText(rv$step6_stats_text)
  output$step6_table <- renderDT({
    req(rv$filter_results)
    datatable(head(rv$filter_results$filtered_data, 200),
              options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step6_fusion_table <- renderDT({
    req(rv$fusion_summary)
    datatable(rv$fusion_summary, options = list(dom = 't'))
  })

  # ===========================================================================
  #  STEP 7  —  Hybrid ML Pipeline
  # ===========================================================================
  observeEvent(input$run_step7, {
    req(rv$gc_clean)
    rv$log_step7 <- ""
    tryCatch({
      option <- input$step7_option
      datasets_to_run <- input$step7_datasets

      # Validate selections
      if (length(input$step7_impute) == 0 || length(input$step7_norm) == 0 ||
          length(input$step7_algo) == 0) {
        append_log("log_step7", "✗ Please select at least one imputation, normalization, and algorithm.")
        return()
      }
      if (length(datasets_to_run) == 0) {
        append_log("log_step7", "✗ Please select at least one dataset.")
        return()
      }

      append_log("log_step7", sprintf("── Hybrid Pipeline Analysis [%s] ──", option))

      # Build data combinations from cleaned datasets
      all_clean <- list(
        gc      = rv$gc_clean,
        gc_icp  = rv$gc_icp_clean
      )
      data_combinations <- all_clean[names(all_clean) %in% datasets_to_run]

      if (length(data_combinations) == 0) {
        append_log("log_step7", "✗ No datasets selected or available.")
        return()
      }

      # Parallel setup
      if (input$step7_parallel) {
        num_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
        cl <- parallel::makePSOCKcluster(num_cores)
        doParallel::registerDoParallel(cl)
        append_log("log_step7", sprintf("  Parallel: %d cores", num_cores))
      }

      t0 <- Sys.time()
      results <- list()

      withProgress(message = "Running hybrid analysis...", value = 0, {
        n_datasets <- length(data_combinations)
        for (i in seq_along(data_combinations)) {
          name <- names(data_combinations)[i]
          incProgress(1 / n_datasets, detail = paste("Dataset:", name))
          append_log("log_step7", sprintf("\n  Processing: %s ...", name))

        # Determine the split argument based on option
        use_sb_only   <- (option == "sb_only")
        use_source_sp <- (option == "sb_env")

        args <- list(
          data = data_combinations[[name]],
          type_col = "Plastic_type",
          remove_cols = c("Source", "Polymer", "technique"),
          data_name = name,
          imputation_methods = input$step7_impute,
          normalization_methods = input$step7_norm,
          algorithms = input$step7_algo,
          cv_folds = input$step7_cv,
          top_k = input$step7_topk,
          top_percent = input$step7_top_pct,
          stability_threshold = input$step7_stab,
          mean_rank_percentile = input$step7_rank,
          n_permutations = input$step7_nperm,
          parallel = FALSE,
          seed = input$step7_seed,
          min_sample_number = 2,
          verbose = TRUE
        )

        if (use_sb_only) {
          args$use_store_vs_environmental_split <- TRUE
        } else {
          args$use_source_split <- TRUE
        }

        results[[name]] <- do.call(run_hybrid_analysis_manuscript1, args)
        gc()
        }
      }) # end withProgress

      t1 <- Sys.time()
      append_log("log_step7", sprintf("\n  Total run time: %.2f minutes",
                                      as.numeric(difftime(t1, t0, units = "mins"))))

      if (input$step7_parallel) {
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      }

      # Store results
      if (option == "sb_only") {
        rv$hybrid_results_sb <- results
      } else {
        rv$hybrid_results_sb_env <- results
      }

      append_log("log_step7", "✓ Hybrid analysis complete. Check the results tabs.")
    },
    error = function(e) append_log("log_step7", paste("✗ ERROR:", conditionMessage(e))))
  })

  # Step 7 output renderers
  output$step7_log <- renderText(rv$log_step7)

  get_step7_result <- reactive({
    res <- if (input$step7_option == "sb_only") rv$hybrid_results_sb else rv$hybrid_results_sb_env
    if (is.null(res) || length(res) == 0) return(NULL)
    # Prefer the first selected dataset that has results
    for (nm in input$step7_datasets) {
      if (!is.null(res[[nm]])) return(res[[nm]])
    }
    # Fallback to first result
    res[[1]]
  })

  output$step7_best <- renderText({
    r <- get_step7_result()
    req(r)
    paste(
      sprintf("Best Pipeline:"),
      sprintf("  Imputation:    %s", r$best_imputation %||% "N/A"),
      sprintf("  Normalization: %s", r$best_normalization %||% "N/A"),
      sprintf("  Algorithm:     %s", r$best_algorithm %||% "N/A"),
      "",
      sprintf("MCC:               %.4f", r$mcc %||% NA),
      sprintf("Balanced Accuracy: %.4f", r$metrics$macro_recall %||% NA),
      sprintf("Macro F1:          %.4f", r$metrics$macro_f1 %||% NA),
      sep = "\n"
    )
  })

  output$step7_metrics_table <- renderDT({
    r <- get_step7_result()
    req(r)
    m <- r$metrics
    df <- data.frame(
      Metric = c("MCC", "Accuracy", "Balanced Accuracy", "Macro F1", "Kappa"),
      Value  = round(c(m$mcc %||% NA, m$accuracy %||% NA,
                       m$macro_recall %||% NA, m$macro_f1 %||% NA,
                       m$kappa %||% NA), 4)
    )
    datatable(df, options = list(dom = 't'))
  })

  output$step7_cm_plot <- renderPlot({
    r <- get_step7_result()
    req(r, r$plots$confusion_matrix)
    r$plots$confusion_matrix
  })

  output$step7_tier_plot <- renderPlot({
    r <- get_step7_result()
    req(r, r$plots$feature_tiers)
    r$plots$feature_tiers
  })

  output$step7_features_table <- renderDT({
    r <- get_step7_result()
    req(r, r$feature_tiers)
    datatable(r$feature_tiers[, c("feature", "stability_score", "mean_rank", "tier")],
              options = list(scrollX = TRUE, pageLength = 20))
  })

  # ===========================================================================
  #  STEP 8  —  HCA
  # ===========================================================================
  observeEvent(input$run_step8, {
    rv$log_step8 <- ""
    tryCatch({
      append_log("log_step8", "── Hierarchical Cluster Analysis ──")

      # Get the right results object
      src <- input$hca_source
      ds  <- input$hca_dataset
      res_list <- if (src == "sb") rv$hybrid_results_sb else rv$hybrid_results_sb_env
      req(res_list)

      target_obj <- res_list[[ds]]
      if (is.null(target_obj)) target_obj <- res_list[[1]]
      if (is.null(target_obj)) stop("No results found from Step 7.")

      df <- dplyr::bind_rows(
        target_obj$final_imp_norm_train,
        target_obj$final_imp_norm_test
      ) %>%
        dplyr::mutate(Plastic_type = gsub("..", " ", Plastic_type, fixed = TRUE)) %>%
        dplyr::mutate(Plastic_type = gsub(".", " ", Plastic_type, fixed = TRUE))

      if ("Source" %in% names(df)) {
        df <- df %>% dplyr::mutate(Plastic_type_Source = paste0(Plastic_type, "-", Source))
      } else {
        df <- df %>% dplyr::mutate(Plastic_type_Source = Plastic_type)
      }

      hc_df <- df %>%
        dplyr::select(-any_of(c("Subcategory","Source","Plastic_type","Polymer","technique"))) %>%
        dplyr::select(where(is.numeric)) %>%
        as.data.frame()

      hc_df[!is.finite(as.matrix(hc_df))] <- 0
      keep <- rowSums(hc_df, na.rm = TRUE) > 0
      hc_df_clean <- hc_df[keep, , drop = FALSE]
      labels_clean <- df$Plastic_type_Source[keep]

      dist_method <- input$hca_dist
      link_method <- input$hca_link

      if (dist_method == "bray") {
        if (any(hc_df_clean < 0, na.rm = TRUE)) {
          append_log("log_step8", "  Negative values detected — falling back to manhattan.")
          dist_method <- "manhattan"
        }
      }

      d <- if (dist_method == "bray") {
        vegan::vegdist(hc_df_clean, method = "bray")
      } else if (dist_method == "manhattan") {
        vegan::vegdist(hc_df_clean, method = "manhattan")
      } else {
        dist(hc_df_clean, method = dist_method)
      }

      if (any(!is.finite(d))) stop("Distance matrix has non-finite entries.")

      hca_obj <- hclust(d, method = link_method)
      ccc <- cor(d, stats::cophenetic(hca_obj))

      rv$hca_result <- list(hca = hca_obj, labels = labels_clean, ccc = ccc)

      append_log("log_step8", sprintf("  Distance: %s  |  Linkage: %s", dist_method, link_method))
      append_log("log_step8", sprintf("  Cophenetic correlation: %.4f", ccc))
      append_log("log_step8", "✓ HCA complete.")
    },
    error = function(e) append_log("log_step8", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step8_log <- renderText(rv$log_step8)
  output$step8_ccc <- renderText({
    req(rv$hca_result)
    sprintf("Cophenetic Correlation Coefficient (CCC): %.4f", rv$hca_result$ccc)
  })
  output$step8_dendro <- renderPlot({
    req(rv$hca_result)
    par(cex = 0.7, mar = c(5, 4, 2, 1))
    plot(rv$hca_result$hca, labels = rv$hca_result$labels, hang = -1,
         main = "Hierarchical Cluster Analysis", xlab = "", sub = "")
  })

  output$download_hca <- downloadHandler(
    filename = function() paste0("PlastiPrint_HCA_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(rv$hca_result)
      pdf(file, width = 16, height = 8)
      par(cex = 0.6)
      plot(rv$hca_result$hca, labels = rv$hca_result$labels, hang = -1,
           main = "PlastiPrint — HCA Dendrogram", xlab = "", sub = "")
      dev.off()
    }
  )

} # end server

# =============================================================================
#  Launch
# =============================================================================
shinyApp(ui = ui, server = server)
