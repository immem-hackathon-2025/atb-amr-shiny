# plot candidate core genes for a selected species and core threshold
# take inputs from core_gene_species_controls.R

coreGenesForGenusPage <- function(afp, input, output) {
  
  ui <- fluidPage(
      # select Genus for core_gene_species plot
      selectInput(
        "selected_genus",
        "Choose a genus to explore gene frequency across its species:",
        list("Enterobacter"="Enterobacter", "Enterobacter hormaechei"="Enterobacter hormaechei")
      ),
      # select core gene threshold for core_gene_species plot
      sliderInput(
        "core_threshold",
        "Select a minimum frequency threshold for core genes:",
        min=0,max=1,value=0.9
      ),
      # select core gene threshold for core_gene_species plot
      sliderInput(
        "min_genomes_per_species",
        "Select a minimum number of genomes per species, to include in plot:",
        min=5,max=100,value=10
      ),
      plotOutput("coreGeneGenusPlot")
  )


  output$coreGeneGenusPlot <- renderPlot({
    
    # for a single species, plot candidate core genes
    # TODO: user settable
  
    afp <- afp %>% filter(grepl(input$selected_genus, Species))
  
    # total number per species
    n_per_species <- afp %>% 
      group_by(Name, Species) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(Species) %>% 
      summarise(nspp=n()) %>% 
      arrange(-nspp)
    
    ### TODO: allow user to select node instead of gene, as the unit of measurement
    # gene frequency per species
    gene_count_per_spp <- afp %>% 
      group_by(Name, `Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(`Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
      count() %>% 
      left_join(n_per_species, by="Species") %>% 
      mutate(freq=n/nspp)
    
    # extract core genes and display freq per species
    gene_count_per_spp %>% filter(nspp>input$min_genomes_per_species & freq>input$core_threshold) %>%
      mutate(label=paste0(Species, " (n=", nspp, ")")) %>%
      arrange(-nspp) %>%
      ggplot(aes(y=label, x=freq)) + 
      geom_col(position='dodge') + 
      facet_wrap(~`Gene symbol`) + 
      theme(legend.position="bottom")
      
  }) 

  return(ui)

}