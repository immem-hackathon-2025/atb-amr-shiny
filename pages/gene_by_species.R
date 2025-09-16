# plot gene frequency across species

geneAcrossSpeciesPage <- function(afp, input, output) {
  
  # total number per gene
  n_per_gene <- afp %>% 
    group_by(Name, `Gene symbol`) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(`Gene symbol`) %>% 
    summarise(n=n()) 
  
  # add total species per gene
  n_per_gene <- afp %>% 
    group_by(`Gene symbol`, Species) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(`Gene symbol`) %>% 
    summarise(nspp=n()) %>% 
    left_join(n_per_gene) %>%
    arrange(-n) 
  
  gene_list <- n_per_gene$`Gene symbol`
  names(gene_list) <- paste(n_per_gene$`Gene symbol`, " (n=",n_per_gene$n," in ",n_per_gene$nspp," species)", sep="")
  
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        selectInput(
          "selected_gene",
          "Choose a gene to explore its frequency across species:",
          gene_list
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
        shiny::uiOutput("geneAcrossSpeciesDownloadButton"),
      )
    )
  )
  
  geneAcrossSpecies <- reactive({
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
      filter(freq>=input$min_freq & n>input$min_genomes_per_species)
  })
  
  output$geneAcrossSpeciesDownloadButton <- renderUI({
    IconButton("downloadGeneAcrossSpecies", "data_dl", "Download")
  })
  output$downloadGeneAcrossSpecies <- downloadHandler(
    filename = "genes_by_species.tsv",
    content = function(file) {
      write_tsv(geneAcrossSpecies(), file)
    }
  )
  
  output$geneAcrossSpeciesPlot <- renderPlot({
    geneAcrossSpecies() %>% 
      ggplot(aes(x=freq, y=Species)) +
      geom_col(fill="navy") + 
      theme_bw() +
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Species with >=",input$min_genomes_per_species, " genomes and gene ",input$selected_gene," freq >=",input$min_freq))
    
  }) 
  
  return(ui)
  
}