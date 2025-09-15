# plot gene frequency for a selected species and core threshold

coreGenesForSpeciesPage <- function(afp, input, output) {
  
  ui <- fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "selected_species",
            "Choose a species to explore its gene frequency:",
            list("Enterobacter cloacae"="Enterobacter cloacae", "Enterobacter hormaechei"="Enterobacter hormaechei")
          ),
          # select gene threshold for core_gene_species plot
          sliderInput(
            "core_threshold",
            "Select a minimum frequency threshold for genes to include:",
            min=0,max=1,value=0.9
          )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("coreGeneSpeciesPlot"),
        )
      )
  )


  output$coreGeneSpeciesPlot <- renderPlot({
    
    # for a single species, plot candidate core genes
    afp_this_spp <- afp %>% filter(Species==input$selected_species)
    n_this_spp <- length(unique(afp_this_spp$Name))
  
    afp_this_spp %>%
      group_by(Name, `Gene symbol`, Class, Subclass) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(`Gene symbol`, Class, Subclass) %>% 
      count() %>% 
      mutate(freq=n/length(unique(afp_this_spp$Name))) %>%
      filter(freq>=input$core_threshold) %>% 
      ggplot(aes(x=freq, y=`Gene symbol`, fill=Class)) +
      geom_col() + 
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Genes with freq >= ",input$core_threshold," in ", input$selected_species))
      
  }) 
  
  return(ui)

}