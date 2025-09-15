# plot candidate core genes for a selected species and core threshold
# take inputs from core_gene_species_controls.R

coreGenesForSpeciesPage <- function(afp, input, output) {
  
  ui <- fluidPage(
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
      plotOutput("coreGeneSpeciesPlot")
  )


  output$coreGeneSpeciesPlot <- renderPlot({
    
    # for a single species, plot candidate core genes
    # TODO: user settable
  
    afp_this_spp <- afp %>% filter(Species==input$selected_species)
    n_this_spp <- length(unique(afp_this_spp$Name))
    
    # TODO: make pretty, annotate with class/subclass
    afp_this_spp %>%
      group_by(Name, `Gene symbol`, Class, Subclass) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(`Gene symbol`, Class, Subclass) %>% 
      count() %>% 
      mutate(freq=n/length(unique(afp_this_spp$Name))) %>%
      filter(freq>=input$core_threshold) %>% 
      ggplot(aes(x=freq, y=`Gene symbol`)) +
      geom_col()
      
  }) 

  return(ui)

}