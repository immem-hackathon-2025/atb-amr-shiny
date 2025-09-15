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


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

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
    
      plotCoreGenesForSpecies(core_threshold=input$core_threshold, selected_species=input$selected_species)
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
