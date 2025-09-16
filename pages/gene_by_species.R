# plot gene frequency across species

geneAcrossSpeciesPage <- function(afp, input, output) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        selectInput(
          "selected_gene",
          "Choose a gene to explore its frequency across species:",
          list("oqxA"="oqxA", "qnrB1"="qnrB1")
        ),
        
        sliderInput(
          "min_genomes_per_species",
          "Select a minimum number of genomes per species, to include the species in the plot:",
          min=5,max=100,value=10
        ),
        
        sliderInput(
          "min_freq",
          "Select a minimum gene frequency per species, to include the species in the plot:",
          min=0,max=1,value=0.01
        )
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("geneAcrossSpeciesPlot", height = "calc(100vh - 200px)"),
      )
    )
  )
  
  
  output$geneAcrossSpeciesPlot <- renderPlot({
    
    # for a single species, plot candidate core genes
    afp_this_gene <- afp %>% filter(`Gene symbol`==input$selected_gene)
    #n_this_spp <- length(unique(afp_this_spp$Name))
    
    # total number per species
    n_per_species <- afp %>% 
      group_by(Name, Species) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(Species) %>% 
      summarise(nspp=n()) %>% 
      arrange(-nspp)
    
    afp_this_gene %>%
      group_by(Name, Species, Class, Subclass) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(Species, Class, Subclass) %>% 
      count() %>% 
      left_join(n_per_species, by="Species") %>% 
      mutate(freq=n/nspp) %>%
      filter(freq>=input$min_freq & n>input$min_genomes_per_species) %>% 
      ggplot(aes(x=freq, y=Species)) +
      geom_col(fill="navy") + 
      theme_bw() +
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Species with >=",input$min_genomes_per_species, " genomes and gene ",input$selected_gene," freq >=",input$min_freq))
    
  }) 
  
  return(ui)
  
}