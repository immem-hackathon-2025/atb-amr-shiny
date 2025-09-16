# plot gene frequencies for a selected species and threshold

coreGenesForGenusPage <- function(afp, input, output) {
  
  afp <- afp %>% 
    separate(col = Species, into = c("Genus", "species_name"), sep = " ", remove=F)
  
  # get list of genera to select from
  # total number per genus
  n_per_genus <- afp %>% 
    group_by(Name, Genus) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(Genus) %>% 
    count() %>% 
    arrange(-n)
  
  # add total species per genus
  n_per_genus <- afp %>% 
    group_by(Genus, Species) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(Genus) %>% 
    summarise(nspp=n()) %>% 
    left_join(n_per_genus, by = "Genus") %>%
    arrange(-nspp) 
  
  genus_list <- n_per_genus$Genus
  names(genus_list) <- paste(n_per_genus$Genus, " (n=",n_per_genus$n," in ", n_per_genus$nspp," species)", sep="")
  
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
            shiny::uiOutput("geneCountPerSppDownloadButton"),
          )
      )
  )
  
  geneCountPerSpp <- reactive({
    # for a single species, plot candidate core genes
    #afp <- afp %>% filter(grepl(input$selected_genus, Species))
    afp_this_genus <- afp %>% filter(Genus==input$selected_genus)
    
    # total number per species
    n_per_species <- afp_this_genus %>% 
      group_by(Name, Species) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(Species) %>% 
      summarise(nspp=n()) %>% 
      arrange(-nspp)
    
    ### TODO: allow user to select node instead of gene, as the unit of measurement
    # gene frequency per species
    afp_this_genus %>% 
      group_by(Name, `Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
      count() %>% distinct() %>% ungroup() %>% 
      group_by(`Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
      count() %>% 
      left_join(n_per_species, by="Species") %>% 
      mutate(freq=n/nspp) %>% 
      filter(nspp>input$min_genomes_per_species & freq>input$core_threshold2) %>%
      mutate(label=paste0(Species, " (n=", nspp, ")")) %>%
      arrange(-nspp) 
  })
  
  output$geneCountPerSppDownloadButton <- renderUI({
    IconButton("downloadGeneCountPerSpp", "data_dl", "Download")
  })
  output$downloadGeneCountPerSpp <- downloadHandler(
    #filename = "gene_count_per_spp.tsv",
    filename = paste0("gene_count_per_", input$selected_genus, "_spp.tsv", sep=""),
    content = function(file) {
      write_tsv(geneCountPerSpp(), file)
    }
  )

  output$coreGeneGenusPlot <- renderPlot({
    
    if (nrow(geneCountPerSpp()) > 0) {
      geneCountPerSpp() %>%
      ggplot(aes(y=label, x=freq)) + 
      geom_col(position='dodge', fill="navy") + 
      facet_wrap(~`Gene symbol`) + 
      theme_bw() +
      theme(legend.position="bottom")
    }
    
    else {
      ggplot() +
        # Print a message to the plot area
        annotate(
          "text", 
          x = 0.5, 
          y = 0.5, 
          label = "No data passess current filters, try again",
          size = 6, 
          color = "gray40"
        ) +
        # Remove axes and labels to make the plot completely blank
        theme_void()
    }
      
  }) 
  
  

  return(ui)

}