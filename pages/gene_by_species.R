# plot gene frequency across species

geneAcrossSpeciesPage <- function(afp, input, output) {
  
  gene_counts_lazy <- afp %>%
    group_by(`Gene symbol`) %>%
    summarise(
      n = n(),
      nspp = n_distinct(Species)
    ) %>%
    arrange(desc(n))
  
  n_per_gene <- gene_counts_lazy %>%
    collect()
  
  gene_list <- n_per_gene$`Gene symbol`
  names(gene_list) <- paste(
    n_per_gene$`Gene symbol`,
    " (n=", n_per_gene$n,
    " in ", n_per_gene$nspp,
    " species)",
    sep = ""
  )
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        selectInput(
          "selected_gene",
          "Choose a gene to explore its frequency across species:",
          gene_list
        ),
        # Settings button
        dropdownButton(
          tags$h3("Settings"),
        
        sliderInput(
          "min_genomes_per_species",
          "Select a minimum number of genomes per species, to include the species in the plot:",
          min=5,max=100,value=10
        ),
        
        sliderInput(
          "min_freq",
          "Select a minimum gene frequency per species, to include the species in the plot:",
          min=0,max=1,value=0.01
        ),
        #setting button options
        circle = TRUE,
        status = "warning", 
        icon = icon("gear"), width = "300px",
        tooltip = tooltipOptions(title = "Click to change settings"))
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
    
    req(input$selected_gene, input$min_genomes_per_species, input$min_freq)
    
    afp_this_gene <- afp %>% filter(`Gene symbol` == !!input$selected_gene)
    
    # total number per species
    species_counts <- afp %>%
      distinct(Name, Species) %>%                  # unique sample-species pairs
      group_by(Species) %>%
      summarise(nspp = n(), .groups = "drop")
    
    print(species_counts)
    
    afp_this_gene <- afp_this_gene %>%
      distinct(Name, Species, Class, Subclass) %>%
      group_by(Species, Class, Subclass) %>%
      count() %>%
      left_join(species_counts, by="Species") %>%
      mutate(freq=n/nspp) %>%
      filter(freq>=!!input$min_freq) %>%
      filter(n>!!input$min_genomes_per_species) %>%
      collect()
    
    print(afp_this_gene)
    
    afp_this_gene
    
  })
  
  output$geneAcrossSpeciesDownloadButton <- renderUI({
    IconButton("downloadGeneAcrossSpecies", "data_dl", "Download")
  })
  output$downloadGeneAcrossSpecies <- downloadHandler(
    filename = paste0(input$selected_gene, "_distribution.tsv", sep=""),
    content = function(file) {
      write_tsv(geneAcrossSpecies(), file)
    }
  )
  
  output$geneAcrossSpeciesPlot <- renderPlot({
    
    df <- geneAcrossSpecies()
    
    validate(need(nrow(df) > 0, "No species with this gene pass current thresholds."))
    
    ggplot(df, aes(x=freq, y=Species)) +
      geom_col(fill="navy") + 
      theme_bw() +
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Species with >=",input$min_genomes_per_species, " genomes and gene ",input$selected_gene," freq >=",input$min_freq))
    
  }) 
  
  return(ui)
  
}