# plot gene frequency for a selected species and core threshold

coreGenesForSpeciesPage <- function(afp, input, output) {
  
  # total number per species
  n_per_species <- afp %>% 
    group_by(Name, Species) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(Species) %>% 
    summarise(nspp=n()) %>% 
    arrange(-nspp) %>% 
    filter(nspp>=10) # hard code to only show species with at least n=10 members
  
  species_list <- n_per_species$Species
  names(species_list) <- paste(n_per_species$Species, " (n=",n_per_species$nspp,")", sep="")

  ui <- fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "selected_species",
            "Choose a species to explore its gene frequency:",
            species_list
            #list("Enterobacter cloacae"="Enterobacter cloacae", "Enterobacter hormaechei"="Enterobacter hormaechei")
          ),
          # select gene threshold for core_gene_species plot
          sliderInput(
            "core_threshold",
            "Select a minimum frequency threshold for genes to include:",
            min=0,max=1,value=c(0.6, 1)
          )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          shiny::fluidRow( 
            align = "center",
            shiny::column(
              width = 12,
              h3("Core genes by species"),
              uiOutput("coreGeneSpeciesDownloadButton")
            ),
          ),
          plotOutput("coreGeneSpeciesPlot", height = "calc(100vh - 200px)"),
        )
      )
  )

  coreGeneSpecies <- reactive({
    afp_this_spp <- afp %>% filter(Species==input$selected_species)
    n_this_spp <- length(unique(afp_this_spp$Name))
    afp_this_spp %>%
      group_by(Name, `Gene symbol`, Class, Subclass) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(`Gene symbol`, Class, Subclass) %>% 
      count() %>% 
      mutate(freq=n/length(unique(afp_this_spp$Name))) %>%
      filter(freq >= input$core_threshold[1] & freq <= input$core_threshold[2])
  })
  
  output$coreGeneSpeciesDownloadButton <- renderUI({
    IconButton("downloadCoreGeneSpecies", "data_dl")
  })
  output$downloadCoreGeneSpecies <- downloadHandler(
    filename = "core_gene_species.tsv",
    content = function(file) {
      d <- coreGeneSpecies()
      write_tsv(coreGeneSpecies(), file)
    }
  )

  output$coreGeneSpeciesPlot <- renderPlot({
    coreGeneSpecies() %>% 
      ggplot(aes(x=freq, y=`Gene symbol`, fill=Class)) +
      geom_col() + 
      theme_bw() +
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Genes with frequency in the range [",input$core_threshold[1],"-", input$core_threshold[2],"] in ", input$selected_species))
      
  }) 
  
  return(ui)

}