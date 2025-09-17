# plot gene frequency across species

geneAcrossSpeciesPage <- function(afp, input, output) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        # filter for AMR/STRESS/VIRULENCE
        shinyWidgets::checkboxGroupButtons(
          inputId = "element_type",
          label = "Determinant(s) of interest",
          choices = c("AMR", "VIRULENCE", "STRESS"), 
          selected= c("AMR")
          ),
        
        selectizeInput(
          "selected_gene",
          "Choose a gene to explore its frequency across species:",
          choices = NULL,
          options = list(
            placeholder = 'Select a gene...',
            onInitialize = I('function() { this.setValue(""); }')
          )
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
        sliderInput(
          "identity_threshold",
          "Minimum nucleotide identity:",
          min = 0.5, max = 1.0, value = 0.9
        ),
        sliderInput(
          "coverage_threshold",
          "Minimum coverage:",
          min = 0.5, max = 1.0, value = 0.9
        ),
        checkboxInput(
          "exclude_partial", 
          "Exclude partial hits", 
          value = TRUE
        ),
        
        #setting button options
        circle = TRUE,
        status = "warning", 
        icon = icon("gear"), width = "300px",
        tooltip = tooltipOptions(title = "Click to change settings")
        )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("geneAcrossSpeciesPlot", height = "calc(100vh - 200px)"),
        shiny::uiOutput("geneAcrossSpeciesDownloadButton"),
      )
    )
  )
  
  filtered_gene_list <- reactive({
    req(input$element_type)
    
    filtered <- afp %>%
      filter(`Element%20type` %in% input$element_type) %>%
      group_by(`Gene symbol`) %>%
      summarise(
        n = n(),
        nspp = n_distinct(Species),
        .groups = "drop"
      ) %>%
      filter(!is.na(`Gene symbol`)) %>%
      arrange(desc(n)) %>%
      collect()
    
    genes <- filtered$`Gene symbol`
    names(genes) <- paste0(
      filtered$`Gene symbol`,
      " (n=", filtered$n,
      " in ", filtered$nspp,
      " species)"
    )
    
    genes
  })
  
  # Server-side selectize for gene choices
  observe({
    gene_choices <- filtered_gene_list()
    
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "selected_gene",
      choices = gene_choices,
      selected = if (length(gene_choices)) gene_choices[[1]] else NULL,
      server = TRUE
    )
  })

  geneAcrossSpecies <- reactive({
    
    
    # for a single species, plot candidate core genes
    
    req(input$selected_gene, input$min_genomes_per_species, input$min_freq, input$identity_threshold, input$coverage_threshold)
    
    # total number per species
    species_counts <- afp %>%
      distinct(Name, Species) %>%                  # unique sample-species pairs
      group_by(Species) %>%
      summarise(nspp = n(), .groups = "drop")
    
    # Convert proportions [0,1] to percentages [0,100] to match the AFP columns
    id_min  <- input$identity_threshold  * 100
    cov_min <- input$coverage_threshold  * 100
    
    afp_this_gene <- afp %>% 
      filter(`Element%20type` %in% !!input$element_type) %>%
      filter(`Gene symbol` == !!input$selected_gene) %>%
      filter(`% Coverage of reference sequence` >= !!cov_min) %>%
      filter(`% Identity to reference sequence` >= !!id_min)
    
    if (input$exclude_partial) {
      afp_this_gene <- afp_this_gene %>% filter(!grepl("PARTIAL", Method))
    }
    
    afp_this_gene <- afp_this_gene %>%
      distinct(Name, Species, Class, Subclass) %>%
      group_by(Species, Class, Subclass) %>%
      count() %>%
      left_join(species_counts, by="Species") %>%
      mutate(freq=n/nspp) %>%
      filter(freq>=!!input$min_freq) %>%
      filter(n>!!input$min_genomes_per_species) %>%
      mutate(label=paste0(Species, " (n=", nspp, ")")) %>%
      collect()
    
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
    
    ggplot(df, aes(x=freq, y=label)) +
      geom_col(fill="navy") + 
      theme_bw() +
      theme(axis.text.y=element_text(size=10)) + 
      labs(y="", x="Gene frequency", title=paste0("Species with >=",input$min_genomes_per_species, " genomes and gene ",input$selected_gene," freq >=",input$min_freq))
    
  }) 
  
  return(ui)
  
}