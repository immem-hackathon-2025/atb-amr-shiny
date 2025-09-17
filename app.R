#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(tidyverse)
library(dplyr)
library(bslib)
library(DBI)
library(duckdb)
library(ggplot2)
library(readr)
library(shinyWidgets)

source("pages/core_gene_species.R")
source("pages/core_gene_genus.R")
source("pages/gene_by_species.R")

source("src/download_functions.R")


# read in raw data, this has a merge of AFP and the species call + HQ
# Connect (single per-page)
con <- dbConnect(duckdb())

# Point to your Parquet dataset (Hive-style partitioned by Genus)
afp_amr <- tbl(con, "read_parquet('data/amr_by_genus/**/*.parquet', hive_partitioning = true)")
afp_stress <- tbl(con, "read_parquet('data/stress_by_genus/**/*.parquet', hive_partitioning = true)")
afp_virulence <- tbl(con, "read_parquet('data/virulence_by_genus/**/*.parquet', hive_partitioning = true)")

connections <- list(
  "AMR" = afp_amr,
  "STRESS" = afp_stress,
  "VIRULENCE" = afp_virulence
)

connections <- list(
  "AMR" = afp_amr,
  "STRESS" = afp_stress,
  "VIRULENCE" = afp_virulence
)

# Load precomputed data from cache files
cache_dir <- "data"
message("Loading precomputed data from cache...")

# Check if cache files exist
cache_files <- c("species_data.rds", "genus_data.rds", "gene_data.rds")
missing_files <- !file.exists(file.path(cache_dir, cache_files))

if (any(missing_files)) {
  stop(paste(
    "Cache files missing! Please run precompute_data.R first to generate:",
    paste(cache_files[missing_files], collapse = ", "),
    "\nRun: Rscript precompute_data.R"
  ))
}

# Load cached data
species_data <- readRDS(file.path(cache_dir, "species_data.rds"))
genus_data <- readRDS(file.path(cache_dir, "genus_data.rds"))
gene_data <- readRDS(file.path(cache_dir, "gene_data.rds"))

# Load and display metadata
if (file.exists(file.path(cache_dir, "metadata.rds"))) {
  metadata <- readRDS(file.path(cache_dir, "metadata.rds"))
  message(paste("Cache data computed at:", metadata$computed_at))
  message("Data summary:")
  for (et in metadata$element_types) {
    message(paste0("  ", et, ": ", metadata$species_counts[et], " species, ", 
                   metadata$genus_counts[et], " genera, ", 
                   metadata$gene_counts[et], " genes"))
  }
} else {
  message("Cache metadata not found, but data loaded successfully")
}

message("Precomputed data loaded successfully!")

# Define UI for application that draws a histogram
ui <- page_navbar(
  title = "AllTheBacteria AMR+ Explorer",
  theme = bs_theme(preset = "lumen"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$style(HTML("
        html, body {
          height: 100vh !important;
        }
      "))
  ),
  nav_panel("Species view", uiOutput("core_genes_species"), icon = icon("bacterium")),
  nav_panel("Gene view", uiOutput("gene_distribution"), icon = icon("dna")),
  nav_panel("Common genes by genus", uiOutput("core_genes_genus"), icon = icon("project-diagram"))
  
)

# Define server logic required to draw a core gene plot for a selected species
server <- function(input, output) {
  
  
  # Render pages
  output$core_genes_species <- renderUI({
    coreGenesForSpeciesPage(connections, species_data, input, output)
  })
  output$core_genes_genus <- renderUI({
    coreGenesForGenusPage(connections, genus_data, input, output)
  })
  output$gene_distribution <- renderUI({
    geneAcrossSpeciesPage(connections, gene_data, input, output)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)