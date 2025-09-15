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


source("core_gene_species.R")
source("core_gene_genus.R")


# read in raw data, this has a merge of AFP and the species call + HQ
### TODO: REPLACE WITH INTERNAL DATA OBJECT
### IDEA: could remove filter to AMR, so we can plot the same info for virulence genes etc reported by AMRfp
afp <-read_tsv("ATB_Enterobacter_AFP.tsv.gz") %>% filter(HQ) %>% filter(`Element subtype`=="AMR")


# Define UI for application that draws a histogram
ui <- page_navbar(
  title = "AllTheBacteria AMR Explorer",
  theme = bs_theme(preset = "lumen"),
  tags$head(
    tags$style(HTML("
        html, body {
          height: 100vh !important;
        }
      "))
  ),
  nav_panel("Core Genes", uiOutput("core_genes")),
)

# Define server logic required to draw a core gene plot for a selected species
server <- function(input, output) {
    # Render pages
    output$core_genes <- renderUI({
        coreGenesForSpeciesPage(afp, input, output)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
