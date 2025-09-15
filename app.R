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


source("pages/core_gene_species.R")
source("pages/core_gene_genus.R")
#source("pages/gene_distribution.R")


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
  nav_panel("Genes per Species", uiOutput("core_genes_species"), icon = icon("dna")),
  nav_panel("Genes per Genus", uiOutput("core_genes_genus"), icon = icon("project-diagram"))

)

# Define server logic required to draw a core gene plot for a selected species
server <- function(input, output) {
    # Render pages
    output$core_genes_species <- renderUI({
        coreGenesForSpeciesPage(afp, input, output)
    })
    output$core_genes_genus <- renderUI({
        coreGenesForGenusPage(afp, input, output)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
