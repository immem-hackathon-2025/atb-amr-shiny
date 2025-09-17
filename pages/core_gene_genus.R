# plot gene frequencies for a selected species and threshold

coreGenesForGenusPage <- function(connections, genus_data, input, output) {
  
  # Reactive to get current database connection based on element_type
  current_afp <- reactive({
    req(input$element_type)
    connections[[input$element_type]]
  })
  
  afp_with_genus <- reactive({
    current_afp() %>%
      mutate(
        Genus = sql("substr(Species, 1, instr(Species, ' ') - 1)"),
        species_name = sql("substr(Species, instr(Species, ' ') + 1)")
      )
  })
  
  # Use precomputed genus data instead of calculating on demand
  genus_counts <- reactive({
    req(input$element_type)
    genus_data[[input$element_type]]$counts
  })
  
  genus_list <- reactive({
    req(input$element_type)
    genus_data[[input$element_type]]$choices
  })
  
  ui <- fluidPage(
      sidebarLayout(
          sidebarPanel(
            # Select for AMR/VIRULENCE/STRESS
            shinyWidgets::checkboxGroupButtons(
              inputId = "element_type",
              label = "Determinant(s) of interest",
              choices = c("AMR", "VIRULENCE", "STRESS"), 
              selected= c("AMR")
            ),
        # select Genus for core_gene_species plot
        selectizeInput(
          "selected_genus",
          "Choose a genus to explore gene frequencies across its member species:",
          choices = NULL,
          options = list(
            placeholder = 'Select a genus...',
            onInitialize = I('function() { this.setValue(""); }')
          )
        ),
        # Settings button
        dropdownButton(
          tags$h3("Settings"),
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
        circle = TRUE,
        status = "warning", 
        icon = icon("gear"), width = "300px",
        tooltip = tooltipOptions(title = "Click to change settings"))
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
  
  # Server-side selectize for genus choices - update when element_type changes
  observe({
    req(genus_list())
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "selected_genus",
      choices = genus_list(),
      selected = if (length(genus_list())) genus_list()[[1]] else NULL,
      server = TRUE
    )
  })

  genus_tbl <- reactive({
    req(input$selected_genus, input$element_type)
    afp_with_genus() %>% filter(Genus == !!input$selected_genus)
  })
  
  geneCountPerSpp <- reactive({
    # for a single species, plot candidate core genes
    
    req(input$core_threshold2, input$selected_genus, input$min_genomes_per_species, input$identity_threshold, input$coverage_threshold)
    
    # Convert proportions [0,1] to percentages [0,100] to match the AFP columns
    id_min  <- input$identity_threshold  * 100
    cov_min <- input$coverage_threshold  * 100
    
    # total number per species
    species_counts <- afp_with_genus() %>%
      distinct(Name, Species) %>%                  # unique sample-species pairs
      group_by(Species) %>%
      summarise(nspp = n(), .groups = "drop")
    
    # Convert proportions [0,1] to percentages [0,100] to match the AFP columns
    id_min  <- input$identity_threshold  * 100
    cov_min <- input$coverage_threshold  * 100
    
    afp_this_genus <- genus_tbl() %>%
      filter(`% Coverage of reference sequence` >= !!cov_min) %>%
      filter(`% Identity to reference sequence` >= !!id_min)
    
    if (input$exclude_partial) {
      afp_this_genus <- afp_this_genus %>% filter(!grepl("PARTIAL", Method))
    }
    
    ### TODO: allow user to select node instead of gene, as the unit of measurement
    # gene frequency per species
    afp_this_genus <- afp_this_genus %>%
      filter(!is.na(`Gene symbol`)) %>% 
      distinct(Name, `Gene symbol`, Class, Subclass, Species) %>%
      group_by(`Gene symbol`, Class, Subclass, Species) %>%
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
    filename = function(){
      paste0("gene_count_per_", input$selected_genus, "_spp.tsv", sep="")
    },
    content = function(file) {
      write_tsv(geneCountPerSpp(), file)
    }
  )

  output$coreGeneGenusPlot <- renderPlot({
    
    df <- geneCountPerSpp()
    
    validate(need(nrow(df) > 0, "No species with this gene pass current thresholds."))
    
    ggplot(df, aes(y=label, x=freq)) + 
      geom_col(position='dodge', fill="navy") + 
      facet_wrap(~`Gene symbol`) + 
      theme_bw() +
      theme(legend.position="bottom")
      
  }) 
  
  

  return(ui)

}