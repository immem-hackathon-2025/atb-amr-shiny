# plot gene frequencies for a selected species and threshold

coreGenesForGenusPage <- function(afp, input, output) {
  
  afp <- afp %>%
    mutate(
      Genus = sql("substr(Species, 1, instr(Species, ' ') - 1)"),
      species_name = sql("substr(Species, instr(Species, ' ') + 1)")
    )
  
  # get list of genera to select from
  # total number per genus
  genus_counts <- afp %>%
    distinct(Name, Genus) %>%                  # unique sample-species pairs
    group_by(Genus) %>%
    summarise(n = n(), .groups = "drop") 
  
  genus_counts <- afp %>%
    distinct(Species, Genus) %>%                  # unique sample-species pairs
    group_by(Genus) %>%
    summarise(nspp = n(), .groups = "drop") %>%
    filter(nspp >= 10) %>%
    left_join(genus_counts, by="Genus") %>%
    arrange(desc(nspp)) %>% collect()
  
  genus_list <- genus_counts$Genus
  names(genus_list) <- paste(
    genus_counts$Genus,
    " (n=", genus_counts$n,
    " in ", genus_counts$nspp,
    " species)",
    sep = ""
  )
  
  ui <- fluidPage(
      sidebarLayout(
          sidebarPanel(
        # select Genus for core_gene_species plot
        selectInput(
          "selected_genus",
          "Choose a genus to explore gene frequencies across its member species:",
          genus_list
        ),
        # select core gene threshold for core_gene_species plot
        sliderInput(
          "core_threshold2",
          "Select a minimum gene frequency, to include the gene in the plot:",
          min=0,max=1,value=0.9
        ),
        # select core gene threshold for core_gene_species plot
        sliderInput(
          "min_genomes_per_species",
          "Select a minimum number of genomes per species, to include the species in the plot:",
          min=5,max=100,value=10
        ),
          ),
          mainPanel(
            plotOutput("coreGeneGenusPlot", height = "calc(100vh - 200px)"),
            div(
              class = "centered-items-row",
              div(uiOutput("filteredDataGenusDownloadButton")),
              div(uiOutput("geneCountPerSppDownloadButton")),
            )
          )
      )
  )
  
  genus_tbl <- reactive({
    afp %>% filter(Genus == !!input$selected_genus)
  })
  
  geneCountPerSpp <- reactive({
    # for a single species, plot candidate core genes
    
    req(input$core_threshold2, input$selected_genus, input$min_genomes_per_species)
    
    afp_this_genus <- genus_tbl()
    
    # total number per species
    species_counts <- afp %>%
      distinct(Name, Species) %>%                  # unique sample-species pairs
      group_by(Species) %>%
      summarise(nspp = n(), .groups = "drop")
    
    ### TODO: allow user to select node instead of gene, as the unit of measurement
    # gene frequency per species
    afp_this_genus <- afp_this_genus %>%
      filter(!is.na(`Gene symbol`)) %>% 
      distinct(Name, `Gene symbol`, Class, Subclass, Species, `Element type`) %>%
      group_by(`Gene symbol`, Class, Subclass, Species, `Element type`) %>%
      count() %>% 
      left_join(species_counts, by="Species") %>% 
      mutate(freq=n/nspp) %>% 
      filter(nspp>!!input$min_genomes_per_species, freq>!!input$core_threshold2) %>%
      mutate(label=paste0(Species, " (n=", nspp, ")")) %>%
      arrange(-nspp) %>% 
      collect()
    
    afp_this_genus
  })
  
  
  # Download data
  # Filtered data
  output$filteredDataGenusDownloadButton <- renderUI({
    IconButton("downloadfilteredDataGenus", "data_dl", 
               paste(input$selected_genus, "data"))
  })
  output$downloadfilteredDataGenus <- downloadHandler(
    filename = function() {
      paste0(gsub("\\s+", "_", input$selected_genus), "_AFP_data.tsv")
    },
    content = function(file) {
      write_tsv(genus_tbl() %>% collect(), file)
    }
  )
  # frequencies
  output$geneCountPerSppDownloadButton <- renderUI({
    IconButton("downloadGeneCountPerSpp", "data_dl", "Gene count")
  })
  output$downloadGeneCountPerSpp <- downloadHandler(
    filename = paste0("gene_count_per_", input$selected_genus, ".tsv", sep=""),
    content = function(file) {
      write_tsv(geneCountPerSpp() %>% collect(), file)
    }
  )

  output$coreGeneGenusPlot <- renderPlot({
    
    df <- geneCountPerSpp()
    
    validate(need(nrow(df) > 0, "No species with this gene pass current thresholds."))
    
    #if (nrow(geneCountPerSpp()) > 0) {
      ggplot(df, aes(y=label, x=freq)) + 
      geom_col(position='dodge', fill="navy") + 
      facet_wrap(~`Gene symbol`) + 
      theme_bw() +
      theme(legend.position="bottom")
    #}
    
    # else {
    #   ggplot() +
    #     # Print a message to the plot area
    #     annotate(
    #       "text", 
    #       x = 0.5, 
    #       y = 0.5, 
    #       label = "No data passess current filters, try again",
    #       size = 6, 
    #       color = "gray40"
    #     ) +
    #     # Remove axes and labels to make the plot completely blank
    #     theme_void()
    # }
      
  }) 
  
  

  return(ui)

}