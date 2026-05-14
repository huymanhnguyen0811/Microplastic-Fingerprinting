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
helper_file <- file.path(getwd(), "Helper Function using only RF_Github_13May2026.R")
if (!file.exists(helper_file)) {
  # Also try one level up or in the app directory
  alt <- list.files(path = c(".", ".."), pattern = "Helper.*Function.*\\.R$",
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
      menuItem("Step 6 — Features", tabName = "step6",    icon = icon("cut")) #,
      # menuItem("Step 7 — ML",       tabName = "step7",    icon = icon("brain")),
      # menuItem("Step 8 — HCA",      tabName = "step8",    icon = icon("project-diagram"))
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
              title = "Welcome to PlastiPrint - A R Shiny App for Data processing & Source tracking of Complex contaminant mixtures",
              h4("This GUI was developed from the publication, titled:", 
                 tags$b("Chemical Fingerprints of New vs. Weathered Microplastics")),
              p("This GUI walks you through the complete computational workflow for microplastic chemical fingerprinting. Each tab on the left corresponds to a major step."),
              tags$ol(
                tags$li(tags$b("Step 1:"), " Import raw data (ATD-GC-MS / HPLC-TOF-MS) and group compounds"),
                tags$li(tags$b("Step 2:"), " Blank subtraction & negative-peak removal"),
                tags$li(tags$b("Step 3:"), " Assign source (Store-Bought / Environmental)"),
                tags$li(tags$b("Step 4:"), " Import ICP-MS trace-metal data"),
                tags$li(tags$b("Step 5:"), " Merge label information & create dataset combinations"),
                tags$li(tags$b("Step 6:"), " Feature reduction (Regular 80%, Modified 80%, In-House) & data fusion")
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
              p("Click the button below to select your data files manually."),
              actionButton("open_import_modal", "Import Data",
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

    # Logs
    log_step1 = "Awaiting input...",
    log_step2 = "Run Step 1 first.",
    log_step3 = "Run Step 2 first.",
    log_step4 = "Awaiting input...",
    log_step5 = "Run Steps 1-4 first.",
    log_step6 = "Run Step 5 first.",
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

  # ── 1.1: Show Import Modal Dialog ───────────────────────────────────────────
  observeEvent(input$open_import_modal, {
    showModal(modalDialog(
      title = "Import Data Files",
      size = "l",
      easyClose = FALSE,
      
      fluidRow(
        column(12,
          h4(icon("flask"), " ATD-GC-MS Data (CSV files)"),
          p("Select CSV files for ATD-GC-MS sample data. You can select multiple files."),
          fileInput("gc_sample_files", 
                    label = "GC Sample Files (.csv):",
                    multiple = TRUE,
                    accept = c(".csv", "text/csv"),
                    width = "100%"),
          fileInput("gc_blank_files", 
                    label = "GC Blank Files (.csv) [Optional]:",
                    multiple = TRUE,
                    accept = c(".csv", "text/csv"),
                    width = "100%"),
          hr(),
          
          h4(icon("tint"), " HPLC-TOF-MS Data (XLS files)"),
          p("Select XLS files for HPLC-TOF-MS data. You can select multiple files."),
          fileInput("hplc_files", 
                    label = "HPLC Files (.xls):",
                    multiple = TRUE,
                    accept = c(".xls", ".xlsx", "application/vnd.ms-excel"),
                    width = "100%"),
          hr(),
          
          div(style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px;",
            h5(icon("info-circle"), " File Selection Summary:"),
            textOutput("import_summary")
          )
        )
      ),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_import", "Import Selected Files", 
                     class = "btn-primary", icon = icon("check"))
      )
    ))
  })
  
  # ── Summary of selected files in modal ──────────────────────────────────────
  output$import_summary <- renderText({
    gc_samples <- if (!is.null(input$gc_sample_files)) nrow(input$gc_sample_files) else 0
    gc_blanks <- if (!is.null(input$gc_blank_files)) nrow(input$gc_blank_files) else 0
    hplc <- if (!is.null(input$hplc_files)) nrow(input$hplc_files) else 0
    
    paste0(
      "GC Sample files: ", gc_samples, " | ",
      "GC Blank files: ", gc_blanks, " | ",
      "HPLC files: ", hplc
    )
  })

  # ── 1.1: Import raw data from selected files ─────────────────────────────────
  observeEvent(input$confirm_import, {
    removeModal()
    rv$log_step1 <- ""
    
    tryCatch({
      # --- ATD-GC-MS ---
      append_log("log_step1", "── Importing ATD-GC-MS data ──")
      
      read_and_add_filename <- function(file_path, file_name) {
        data.table::fread(file_path, showProgress = FALSE, data.table = FALSE,
                          check.names = TRUE) %>%
          mutate(File = file_name)
      }
      
      # Process GC Sample Files
      if (!is.null(input$gc_sample_files) && nrow(input$gc_sample_files) > 0) {
        gc_sample_list <- lapply(1:nrow(input$gc_sample_files), function(i) {
          read_and_add_filename(
            input$gc_sample_files$datapath[i],
            gsub(".csv", "", input$gc_sample_files$name[i])
          )
        })
        
        atdgcms_raw <- data.table::rbindlist(gc_sample_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          dplyr::select(-any_of(c("Start","End","Width","Base.Peak",
                                  "Cpd","Label","Height","Ions"))) %>%
          dplyr::mutate(File = gsub("_", "-", File)) %>%
          mutate(type = "Sample")
        
        append_log("log_step1", sprintf("  GC samples: %d rows from %d files",
                                        nrow(atdgcms_raw), nrow(input$gc_sample_files)))
      } else {
        atdgcms_raw <- data.frame()
        append_log("log_step1", "  No GC sample files selected.")
      }
      
      # Process GC Blank Files
      if (!is.null(input$gc_blank_files) && nrow(input$gc_blank_files) > 0) {
        gc_blank_list <- lapply(1:nrow(input$gc_blank_files), function(i) {
          data.table::fread(input$gc_blank_files$datapath[i], 
                           showProgress = FALSE, check.names = TRUE) %>%
            mutate(File = gsub(".csv", "", input$gc_blank_files$name[i]))
        })
        
        atdgcms_blank <- data.table::rbindlist(gc_blank_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          dplyr::select(any_of(c("Area","File","m.z","RT"))) %>%
          mutate(type = "Blanks")
        
        append_log("log_step1", sprintf("  GC blanks : %d rows from %d files",
                                        nrow(atdgcms_blank), nrow(input$gc_blank_files)))
      } else {
        atdgcms_blank <- data.frame()
        append_log("log_step1", "  No GC blank files selected.")
      }
      
      # Combine GC data
      if (nrow(atdgcms_raw) > 0 || nrow(atdgcms_blank) > 0) {
        atdgcms_step1 <- bind_rows(atdgcms_raw, atdgcms_blank) %>%
          arrange(RT) %>%
          mutate(simplified_file = map_chr(File, ~ {
            parts <- str_split_1(.x, pattern = "-")
            if (length(parts) >= 5) {
              paste(parts[c(4, 5)], collapse = "-")
            } else {
              .x
            }
          }))
        
        rv$atdgcms_raw   <- atdgcms_raw
        rv$atdgcms_blank <- atdgcms_blank
        rv$atdgcms_step1 <- atdgcms_step1
      }

      # --- HPLC-TOF-MS ---
      append_log("log_step1", "\n── Importing HPLC-TOF-MS data ──")
      
      if (!is.null(input$hplc_files) && nrow(input$hplc_files) > 0) {
        read_hplc_file <- function(file_path, file_name) {
          df <- readxl::read_xls(file_path, skip = 1) %>% as.data.frame()
          
          # Handle column names
          if ("m/z" %in% names(df)) {
            df <- df %>% dplyr::rename(m.z = `m/z`)
          }
          
          df <- df %>%
            dplyr::select(any_of(c("m.z", "RT", "Height", "Area"))) %>%
            mutate(
              File = gsub("_", "-", gsub("\\.(xls|xlsx)$", "", file_name)),
              type = ifelse(str_detect(file_name, regex("blank", ignore_case = TRUE)), "Blanks", "Sample"),
              Day = case_when(
                str_detect(file_name, "Day0") ~ "0",
                str_detect(file_name, "Day2") ~ "2",
                str_detect(file_name, "Day15") ~ "15",
                TRUE ~ "34"
              ),
              Replicate = ifelse(str_detect(file_name, "R1"), "R1", "R2"),
              Batch_number = "manual"
            )
          return(df)
        }
        
        hplc_list <- lapply(1:nrow(input$hplc_files), function(i) {
          tryCatch({
            read_hplc_file(input$hplc_files$datapath[i], input$hplc_files$name[i])
          }, error = function(e) {
            warning(paste("Error reading file:", input$hplc_files$name[i], "-", e$message))
            return(data.frame())
          })
        })
        
        all_hplc <- data.table::rbindlist(hplc_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          arrange(RT) %>%
          mutate(simplified_file = map_chr(File, ~ {
            parts <- str_split_1(.x, "-")
            if (length(parts) >= 2) {
              paste(parts[1:min(2, length(parts))], collapse = "-")
            } else {
              .x
            }
          }))
        
        rv$all_hplc <- all_hplc
        append_log("log_step1", sprintf("  HPLC total: %d rows from %d files", 
                                        nrow(all_hplc), nrow(input$hplc_files)))
      } else {
        append_log("log_step1", "  No HPLC files selected.")
      }
      
      append_log("log_step1", "\n✓ Data import complete.")
    },
    error = function(e) {
      append_log("log_step1", paste("✗ ERROR:", conditionMessage(e)))
    })
  })

  # ── 1.2: Compound grouping ─────────────────────────────────────────────────
  observeEvent(input$run_step1_group, {
    # Check if at least one dataset is available
    has_gc <- !is.null(rv$atdgcms_step1) && nrow(rv$atdgcms_step1) > 0
    has_hplc <- !is.null(rv$all_hplc) && nrow(rv$all_hplc) > 0
    
    if (!has_gc && !has_hplc) {
      append_log("log_step1", "\n⚠ No data imported yet. Please import data first.")
      return()
    }
    
    tryCatch({
      append_log("log_step1", "\n── Grouping compounds ──")
      
      # Process GC data if available
      if (has_gc) {
        rv$atdgcms_grouped <- grouping_comp(
          rv$atdgcms_step1,
          rtthres1 = input$gc_rtthres1,
          mzthres  = input$gc_mzthres,
          type     = "ATDGCMS")
        append_log("log_step1", sprintf("  GC grouped: %d unique features",
                                        n_distinct(rv$atdgcms_grouped$Feature, na.rm = TRUE)))
      } else {
        append_log("log_step1", "  GC: No data to group (skipped)")
      }

      # Process HPLC data if available
      if (has_hplc) {
        rv$hplc_grouped <- grouping_comp(
          rv$all_hplc,
          split_time = input$hplc_split_time,
          rtthres1   = input$hplc_rtthres1,
          rtthres2   = input$hplc_rtthres2,
          mzthres    = input$hplc_mzthres,
          type       = "HPLCTOFMS")
        append_log("log_step1", sprintf("  HPLC grouped: %d unique features",
                                        n_distinct(rv$hplc_grouped$Feature, na.rm = TRUE)))
      } else {
        append_log("log_step1", "  HPLC: No data to group (skipped)")
      }
      
      append_log("log_step1", "✓ Grouping complete.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  # ── 1.3: Benchmark removal ─────────────────────────────────────────────────
  observeEvent(input$run_step1_bench, {
    # Check if at least one grouped dataset is available
    has_gc <- !is.null(rv$atdgcms_grouped) && nrow(rv$atdgcms_grouped) > 0
    has_hplc <- !is.null(rv$hplc_grouped) && nrow(rv$hplc_grouped) > 0
    
    if (!has_gc && !has_hplc) {
      append_log("log_step1", "\n⚠ No grouped data available. Please run 'Group Compounds' first.")
      return()
    }
    
    tryCatch({
      append_log("log_step1", "\n── Removing benchmarks ──")

      # Process GC data if available
      if (has_gc) {
        gc_g <- rv$atdgcms_grouped
        n_before_gc <- nrow(gc_g)
        idx <- which(gc_g$m.z >= 98 & gc_g$m.z <= 98.5 &
                     gc_g$RT > 13.4 & gc_g$RT < 13.5)
        if (length(idx) > 0) gc_g <- gc_g[-idx, ]
        rv$atdgcms_grouped <- gc_g
        append_log("log_step1", sprintf("  GC: removed %d benchmark rows", n_before_gc - nrow(gc_g)))
      } else {
        append_log("log_step1", "  GC: No data to process (skipped)")
      }

      # Process HPLC data if available
      if (has_hplc) {
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
      } else {
        append_log("log_step1", "  HPLC: No data to process (skipped)")
      }
      
      append_log("log_step1", "✓ Benchmarks removed.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  # Step 1 outputs
  output$step1_log <- renderText(rv$log_step1)
  output$step1_gc_table <- renderDT({
    if (is.null(rv$atdgcms_grouped) || nrow(rv$atdgcms_grouped) == 0) {
      return(datatable(data.frame(Message = "No GC data available. Import and group data first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$atdgcms_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step1_hplc_table <- renderDT({
    if (is.null(rv$hplc_grouped) || nrow(rv$hplc_grouped) == 0) {
      return(datatable(data.frame(Message = "No HPLC data available. Import and group data first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$hplc_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 2  —  Blank Subtraction
  # ===========================================================================
  observeEvent(input$run_step2, {
    # Check if at least one grouped dataset is available
    has_gc <- !is.null(rv$atdgcms_grouped) && nrow(rv$atdgcms_grouped) > 0
    has_hplc <- !is.null(rv$hplc_grouped) && nrow(rv$hplc_grouped) > 0
    
    if (!has_gc && !has_hplc) {
      rv$log_step2 <- "⚠ No grouped data available. Please complete Step 1 first."
      return()
    }
    
    rv$log_step2 <- ""
    tryCatch({
      sd_mult <- input$blank_sd_mult
      append_log("log_step2", sprintf("── Blank subtraction (SD multiplier = %.1f) ──", sd_mult))

      # GC
      if (has_gc) {
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
      } else {
        append_log("log_step2", "  GC: No data to process (skipped)")
      }

      # HPLC
      if (has_hplc) {
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
      } else {
        append_log("log_step2", "  HPLC: No data to process (skipped)")
      }
      
      append_log("log_step2", "✓ Blank subtraction complete.")
    },
    error = function(e) append_log("log_step2", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step2_log <- renderText(rv$log_step2)
  output$step2_gc_table <- renderDT({
    if (is.null(rv$atdgcms_step2) || nrow(rv$atdgcms_step2) == 0) {
      return(datatable(data.frame(Message = "No GC data available. Run blank subtraction first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$atdgcms_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step2_hplc_table <- renderDT({
    if (is.null(rv$hplc_step2) || nrow(rv$hplc_step2) == 0) {
      return(datatable(data.frame(Message = "No HPLC data available. Run blank subtraction first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$hplc_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 3  —  Source assignment
  # ===========================================================================
  observeEvent(input$run_step3, {
    # Check if at least one dataset is available from Step 2
    has_gc <- !is.null(rv$atdgcms_step2) && nrow(rv$atdgcms_step2) > 0
    has_hplc <- !is.null(rv$hplc_step2) && nrow(rv$hplc_step2) > 0
    
    if (!has_gc && !has_hplc) {
      rv$log_step3 <- "⚠ No data available from Step 2. Please complete Step 2 first."
      return()
    }
    
    rv$log_step3 <- ""
    tryCatch({
      append_log("log_step3", "── Assigning Source & Technique ──")

      if (has_gc) {
        gc <- rv$atdgcms_step2 %>%
          dplyr::filter(!is.na(Feature)) %>%
          dplyr::select(File, Feature, m.z, RT, Area) %>%
          mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
          dplyr::rename(Values = Area) %>%
          mutate(technique = "GC")
        rv$gc <- gc
        append_log("log_step3", sprintf("  GC:   %d rows  |  SB: %d  ENV: %d",
          nrow(gc), sum(gc$Source == "Store-Bought"), sum(gc$Source == "Environmental")))
      } else {
        append_log("log_step3", "  GC: No data to process (skipped)")
      }

      if (has_hplc) {
        hplc <- rv$hplc_step2 %>%
          dplyr::filter(!is.na(Feature)) %>%
          dplyr::select(File, Feature, m.z, RT, Height) %>%
          dplyr::mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
          dplyr::rename(Values = Height) %>%
          mutate(technique = "HPLC")
        rv$hplc <- hplc
        append_log("log_step3", sprintf("  HPLC: %d rows  |  SB: %d  ENV: %d",
          nrow(hplc), sum(hplc$Source == "Store-Bought"), sum(hplc$Source == "Environmental")))
      } else {
        append_log("log_step3", "  HPLC: No data to process (skipped)")
      }
      
      append_log("log_step3", "✓ Step 3 complete.")
    },
    error = function(e) append_log("log_step3", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step3_log <- renderText(rv$log_step3)
  output$step3_gc_table <- renderDT({
    if (is.null(rv$gc) || nrow(rv$gc) == 0) {
      return(datatable(data.frame(Message = "No GC data available. Run Step 3 first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$gc, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step3_hplc_table <- renderDT({
    if (is.null(rv$hplc) || nrow(rv$hplc) == 0) {
      return(datatable(data.frame(Message = "No HPLC data available. Run Step 3 first."),
                       options = list(dom = 't')))
    }
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

      icpms_r1 <- read_trace_metal("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/icpms_round1_rawdata_removal USE-01-rep1 and USE-03 (only2observations).xlsx") %>% remove_H2_HMI_modus()
      icpms_r2b1 <- read_trace_metal("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/icpms_round2_batch1_rawdata_26Feb2026.xlsx") %>% remove_H2_HMI_modus()
      icpms_r2b2 <- read_trace_metal("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/icpms_round2_batch2_rawdata.xlsx") %>% remove_H2_HMI_modus()

      rec_unc_1 <- readxl::read_excel("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round1.xlsx")
      rec_unc_2_1 <- readxl::read_excel("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round2_batch1.xlsx")
      rec_unc_2_2 <- readxl::read_excel("./Microplastic manuscript 1 - Raw data/ICPMS_Trace metal/Trace metal data_Recovery_Uncertainty_round2_batch2.xlsx")

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
    # Check what data is available
    has_gc <- !is.null(rv$gc) && nrow(rv$gc) > 0
    has_hplc <- !is.null(rv$hplc) && nrow(rv$hplc) > 0
    has_icp <- !is.null(rv$icp) && nrow(rv$icp) > 0
    
    if (!has_gc && !has_hplc && !has_icp) {
      rv$log_step5 <- "⚠ No data available. Please complete Steps 1-3 (and optionally Step 4 for ICP)."
      return()
    }
    
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

      # Process GC if available
      if (has_gc) {
        gc <- left_join(rv$gc, sampinfo, by = "File")
        rv$gc_labeled <- gc
        append_log("log_step5", sprintf("  GC:   %d rows  (%d plastic types)",
                                        nrow(gc), n_distinct(gc$Plastic_type, na.rm = TRUE)))
      } else {
        gc <- NULL
        append_log("log_step5", "  GC: No data (skipped)")
      }
      
      # Process HPLC if available
      if (has_hplc) {
        hplc <- left_join(rv$hplc, sampinfo, by = "File")
        rv$hplc_labeled <- hplc
        append_log("log_step5", sprintf("  HPLC: %d rows  (%d plastic types)",
                                        nrow(hplc), n_distinct(hplc$Plastic_type, na.rm = TRUE)))
      } else {
        hplc <- NULL
        append_log("log_step5", "  HPLC: No data (skipped)")
      }
      
      # Process ICP if available
      if (has_icp) {
        icp <- left_join(rv$icp, sampinfo, by = "File")
        rv$icp_labeled <- icp
        append_log("log_step5", sprintf("  ICP:  %d rows  (%d plastic types)",
                                        nrow(icp), n_distinct(icp$Plastic_type, na.rm = TRUE)))
      } else {
        icp <- NULL
        append_log("log_step5", "  ICP: No data (skipped)")
      }

      # Create combinations based on available data
      combinations_created <- c()
      
      if (has_gc && has_hplc) {
        gc_hplc <- bind_rows(gc, hplc)
        rv$gc_hplc <- gc_hplc
        combinations_created <- c(combinations_created, "gc_hplc")
      }
      
      if (has_gc && has_icp) {
        gc_icp <- bind_rows(icp, gc %>% dplyr::select(-any_of(c("m.z", "RT"))))
        rv$gc_icp <- gc_icp
        combinations_created <- c(combinations_created, "gc_icp")
      }
      
      if (has_hplc && has_icp) {
        hplc_icp <- bind_rows(icp, hplc %>% dplyr::select(-any_of(c("m.z", "RT"))))
        rv$hplc_icp <- hplc_icp
        combinations_created <- c(combinations_created, "hplc_icp")
      }
      
      if (has_gc && has_hplc && has_icp) {
        gc_hplc_icp <- bind_rows(icp, bind_rows(gc, hplc) %>% dplyr::select(-any_of(c("m.z", "RT"))))
        rv$gc_hplc_icp <- gc_hplc_icp
        combinations_created <- c(combinations_created, "gc_hplc_icp")
      }

      if (length(combinations_created) > 0) {
        append_log("log_step5", sprintf("  Combinations created: %s", paste(combinations_created, collapse = ", ")))
      } else {
        append_log("log_step5", "  No combinations created (need at least 2 data types)")
      }
      
      append_log("log_step5", "✓ Step 5 complete.")
    },
    error = function(e) append_log("log_step5", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step5_log <- renderText(rv$log_step5)
  output$step5_summary_table <- renderDT({
    # Build summary based on available data
    datasets <- c()
    rows <- c()
    features <- c()
    
    if (!is.null(rv$gc_labeled) && nrow(rv$gc_labeled) > 0) {
      datasets <- c(datasets, "gc")
      rows <- c(rows, nrow(rv$gc_labeled))
      features <- c(features, n_distinct(rv$gc_labeled$Feature))
    }
    if (!is.null(rv$hplc_labeled) && nrow(rv$hplc_labeled) > 0) {
      datasets <- c(datasets, "hplc")
      rows <- c(rows, nrow(rv$hplc_labeled))
      features <- c(features, n_distinct(rv$hplc_labeled$Feature))
    }
    if (!is.null(rv$icp_labeled) && nrow(rv$icp_labeled) > 0) {
      datasets <- c(datasets, "icp")
      rows <- c(rows, nrow(rv$icp_labeled))
      features <- c(features, n_distinct(rv$icp_labeled$Feature))
    }
    if (!is.null(rv$gc_hplc) && nrow(rv$gc_hplc) > 0) {
      datasets <- c(datasets, "gc_hplc")
      rows <- c(rows, nrow(rv$gc_hplc))
      features <- c(features, n_distinct(rv$gc_hplc$Feature))
    }
    if (!is.null(rv$gc_icp) && nrow(rv$gc_icp) > 0) {
      datasets <- c(datasets, "gc_icp")
      rows <- c(rows, nrow(rv$gc_icp))
      features <- c(features, n_distinct(rv$gc_icp$Feature))
    }
    if (!is.null(rv$hplc_icp) && nrow(rv$hplc_icp) > 0) {
      datasets <- c(datasets, "hplc_icp")
      rows <- c(rows, nrow(rv$hplc_icp))
      features <- c(features, n_distinct(rv$hplc_icp$Feature))
    }
    if (!is.null(rv$gc_hplc_icp) && nrow(rv$gc_hplc_icp) > 0) {
      datasets <- c(datasets, "gc_hplc_icp")
      rows <- c(rows, nrow(rv$gc_hplc_icp))
      features <- c(features, n_distinct(rv$gc_hplc_icp$Feature))
    }
    
    if (length(datasets) == 0) {
      return(datatable(data.frame(Message = "No data available. Run Step 5 first."),
                       options = list(dom = 't')))
    }
    
    summary_df <- data.frame(Dataset = datasets, Rows = rows, Features = features)
    datatable(summary_df, options = list(dom = 't'))
  })
  output$step5_gc_table <- renderDT({
    if (is.null(rv$gc_labeled) || nrow(rv$gc_labeled) == 0) {
      return(datatable(data.frame(Message = "No GC data available."),
                       options = list(dom = 't')))
    }
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
      
      if (is.null(ds) || nrow(ds) == 0) {
        rv$log_step6 <- sprintf("⚠ Dataset '%s' is not available. Please complete Step 5 first.", ds_name)
        return()
      }
      
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
    # Check what data is available
    has_gc <- !is.null(rv$gc_labeled) && nrow(rv$gc_labeled) > 0
    has_hplc <- !is.null(rv$hplc_labeled) && nrow(rv$hplc_labeled) > 0
    has_icp <- !is.null(rv$icp_labeled) && nrow(rv$icp_labeled) > 0
    
    if (!has_gc && !has_hplc && !has_icp) {
      append_log("log_step6", "\n⚠ No labeled data available. Please complete Step 5 first.")
      return()
    }
    
    tryCatch({
      append_log("log_step6", "\n── Running in-house filter + data fusion ──")
      
      meta <- c("Plastic_type", "technique", "Source", "Polymer")
      fusion_datasets <- c()
      fusion_samples <- c()
      fusion_features <- c()

      # Process individual datasets
      if (has_gc) {
        gc_clean <- shared_feature_filtering(rv$gc_labeled, feat_reduc_mode = "inhouse")$filtered_data
        rv$gc_clean <- gc_clean
        fusion_datasets <- c(fusion_datasets, "gc")
        fusion_samples <- c(fusion_samples, nrow(gc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(gc_clean), meta)))
        append_log("log_step6", sprintf("  GC: %d samples, %d features", nrow(gc_clean), length(setdiff(names(gc_clean), meta))))
      }
      
      if (has_hplc) {
        hplc_clean <- shared_feature_filtering(rv$hplc_labeled, feat_reduc_mode = "inhouse")$filtered_data
        rv$hplc_clean <- hplc_clean
        fusion_datasets <- c(fusion_datasets, "hplc")
        fusion_samples <- c(fusion_samples, nrow(hplc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(hplc_clean), meta)))
        append_log("log_step6", sprintf("  HPLC: %d samples, %d features", nrow(hplc_clean), length(setdiff(names(hplc_clean), meta))))
      }
      
      if (has_icp) {
        icp_clean <- shared_feature_filtering(rv$icp_labeled, feat_reduc_mode = "inhouse")$filtered_data
        rv$icp_clean <- icp_clean
        fusion_datasets <- c(fusion_datasets, "icp")
        fusion_samples <- c(fusion_samples, nrow(icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(icp_clean), meta)))
        append_log("log_step6", sprintf("  ICP: %d samples, %d features", nrow(icp_clean), length(setdiff(names(icp_clean), meta))))
      }

      # Create combinations based on available data
      if (has_gc && has_hplc && !is.null(rv$gc_hplc)) {
        rv$gc_hplc_clean <- data_fusion(rv$gc_hplc, source_feats = list(rv$gc_clean, rv$hplc_clean))
        fusion_datasets <- c(fusion_datasets, "gc_hplc")
        fusion_samples <- c(fusion_samples, nrow(rv$gc_hplc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_hplc_clean), meta)))
      }
      
      if (has_gc && has_icp && !is.null(rv$gc_icp)) {
        rv$gc_icp_clean <- data_fusion(rv$gc_icp, source_feats = list(rv$gc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "gc_icp")
        fusion_samples <- c(fusion_samples, nrow(rv$gc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_icp_clean), meta)))
      }
      
      if (has_hplc && has_icp && !is.null(rv$hplc_icp)) {
        rv$hplc_icp_clean <- data_fusion(rv$hplc_icp, source_feats = list(rv$hplc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "hplc_icp")
        fusion_samples <- c(fusion_samples, nrow(rv$hplc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$hplc_icp_clean), meta)))
      }
      
      if (has_gc && has_hplc && has_icp && !is.null(rv$gc_hplc_icp)) {
        rv$gc_hplc_icp_clean <- data_fusion(rv$gc_hplc_icp, source_feats = list(rv$gc_clean, rv$hplc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "gc_hplc_icp")
        fusion_samples <- c(fusion_samples, nrow(rv$gc_hplc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_hplc_icp_clean), meta)))
      }

      fusion_summary <- data.frame(
        Dataset = fusion_datasets,
        Samples = fusion_samples,
        Features = fusion_features
      )
      rv$fusion_summary <- fusion_summary
      append_log("log_step6", "✓ Data fusion complete. See 'Fusion Summary' tab.")
    },
    error = function(e) append_log("log_step6", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step6_log <- renderText(rv$log_step6)
  output$step6_stats <- renderText(rv$step6_stats_text)
  output$step6_table <- renderDT({
    if (is.null(rv$filter_results)) {
      return(datatable(data.frame(Message = "No filter results available. Run feature reduction first."),
                       options = list(dom = 't')))
    }
    datatable(head(rv$filter_results$filtered_data, 200),
              options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step6_fusion_table <- renderDT({
    if (is.null(rv$fusion_summary)) {
      return(datatable(data.frame(Message = "No fusion summary available. Run data fusion first."),
                       options = list(dom = 't')))
    }
    datatable(rv$fusion_summary, options = list(dom = 't'))
  })
  
  
} # end server

# =============================================================================
#  Launch
# =============================================================================
shinyApp(ui = ui, server = server)
