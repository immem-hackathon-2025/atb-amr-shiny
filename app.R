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

source("core_gene_species.R")


# read in raw data, this has a merge of AFP and the species call + HQ
### TODO: REPLACE WITH INTERNAL DATA OBJECT
### IDEA: could remove filter to AMR, so we can plot the same info for virulence genes etc reported by AMRfp
afp <-read_tsv("ATB_Enterobacter_AFP.tsv.gz") %>% filter(HQ) %>% filter(`Element subtype`=="AMR")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("AllTheBacteria AMR Explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(

        # Show a plot of the generated distribution
        mainPanel(
         
           # select species for core_gene_species plot
          selectInput(
              "selected_species",
              "Choose a species to explore its core genes:",
              list("Enterobacter cloacae"="Enterobacter cloacae", "Enterobacter hormaechei"="Enterobacter hormaechei")
          ),
          
          # select core gene threshold for core_gene_species plot
          sliderInput(
            "core_threshold",
            "Select a minimum frequency threshold for core genes:",
            min=0,max=1,value=0.9
          ),
          
           plotOutput("coreGeneSpecies")
        )
    )
    
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$coreGeneSpecies <- renderPlot({
    
      plotCoreGenesForSpecies(afp=afp, core_threshold=input$core_threshold, selected_species=input$selected_species)
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
