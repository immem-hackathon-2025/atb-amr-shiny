# plot gene frequencies for a selected species and frequency range
coreGenesForSpeciesPage <- function(connections, species_data, input, output) {

  # Reactive to get current database connection based on element_type
  current_afp <- reactive({
    req(input$element_type)
    connections[[input$element_type]]
  })

  # Use precomputed species data instead of calculating on demand
  species_counts <- reactive({
    req(input$element_type)
    species_data[[input$element_type]]$counts
  })
  
  species_choices <- reactive({
    req(input$element_type)
    species_data[[input$element_type]]$choices
  })
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        # Select for AMR/VIRULENCE/STRESS
        shinyWidgets::radioGroupButtons(
          inputId = "element_type",
          label = "Determinant of interest",
          choices = c("AMR", "VIRULENCE", "STRESS"), 
          selected = "AMR"
        ),
        
        selectizeInput(
          "selected_species",
          "Choose a species to explore its gene frequency:",
          choices = NULL,
          options = list(
            placeholder = 'Select a species...',
            onInitialize = I('function() { this.setValue(""); }')
          )
        ),
        
        # Settings button
        dropdownButton(
          tags$h3("Settings"),
          
        # select gene threshold for core_gene_species plot
        sliderInput(
          "core_threshold",
          "Select a frequency range for genes to include:",
          min = 0, max = 1, value = c(0.6, 1.0)
        ),
        # sliders are proportions (0..1); converted below to percent to match columns
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
      mainPanel(
        plotOutput("coreGeneSpeciesPlot", height = "calc(100vh - 200px)"),
        div(
          class = "centered-items-row",
          div(uiOutput("filteredDataSpeciesDownloadButton")),
          div(uiOutput("coreGeneSpeciesDownloadButton")),
        )
      )
    )
  )
  
  # Server-side selectize for species choices - update when element_type changes
  observe({
    req(species_choices())
    updateSelectizeInput(
      session = getDefaultReactiveDomain(),
      inputId = "selected_species",
      choices = species_choices(),
      selected = if (length(species_choices())) species_counts()$Species[[1]] else NULL,
      server = TRUE
    )
  })

  # Reactive data for selected species and thresholds
  spp_tbl <- reactive({
    req(input$selected_species, input$element_type)
    current_afp() %>% filter(Species == !!input$selected_species)
  })
  
  coreGeneSpecies <- reactive({
    req(input$selected_species, input$core_threshold, input$identity_threshold, input$coverage_threshold)
    
    # Convert proportions [0,1] to percentages [0,100] to match the AFP columns
    id_min  <- input$identity_threshold  * 100
    cov_min <- input$coverage_threshold  * 100
    
    # Count how many distinct samples contain each gene (within species)
    # Then compute frequency = n / n_samples_for_species
    # Push everything to DuckDB, collect only the filtered result.
    # Column names with spaces or % need backticks.
    
    n_samples <- spp_tbl() %>%
      summarise(n = n_distinct(Name)) %>%
      collect() %>%
      pull(n)
    
    validate(need(n_samples > 0, paste0("No rows for species: ", input$selected_species)))
    
    gene_freq_tbl <- spp_tbl() %>%
      filter(`% Coverage of reference sequence` >= !!cov_min) %>%
      filter(`% Identity to reference sequence` >= !!id_min)
      
    if (input$exclude_partial) {
      gene_freq_tbl <- gene_freq_tbl %>% filter(!grepl("PARTIAL", Method))
    }
    
    gene_freq_tbl <- gene_freq_tbl %>% 
      filter(!is.na(`Gene symbol`)) %>%
      distinct(Name, `Gene symbol`, Class, Subclass) %>%        # unique sample-gene combos
      group_by(`Gene symbol`, Class, Subclass) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(freq = n / !!n_samples) %>%
      filter(freq >= !!input$core_threshold[1],
             freq <= !!input$core_threshold[2]) %>%
      arrange(desc(freq)) %>%
      collect()
    
    gene_freq_tbl
  })
  
  # Download button
  # Filtered data
  output$filteredDataSpeciesDownloadButton <- renderUI({
    IconButton("downloadfilteredDataSpecies", "data_dl", 
               paste(input$selected_species, "data"))
  })
  output$downloadfilteredDataSpecies <- downloadHandler(
    filename = function() {
      paste0(gsub("\\s+", "_", input$selected_species), "_AFP_data.tsv")
    },
    content = function(file) {
      write_tsv(spp_tbl() %>% collect(), file)
    }
  )
  # Frequencies
  output$coreGeneSpeciesDownloadButton <- renderUI({
    IconButton("downloadCoreGeneSpecies", "data_dl", "Core gene frequencies")
  })
  output$downloadCoreGeneSpecies <- downloadHandler(
    filename = function() {
      paste0("core_genes_", gsub("\\s+", "_", input$selected_species), ".tsv")
    },
    content = function(file) {
      write_tsv(coreGeneSpecies(), file)
    }
  )
  
  # Plot
  output$coreGeneSpeciesPlot <- renderPlot({
    df <- coreGeneSpecies()
    
    validate(need(nrow(df) > 0, "No genes in the selected frequency range with current thresholds."))
    
    basic_plot <- df %>%
      ggplot(aes(x = freq, y = reorder(`Gene symbol`, freq), fill = Class)) +
      geom_col()
    
    # if virulence only, set colour to navy
    if (length(input$element_type)==1) {
      if (input$element_type=="VIRULENCE") {
        basic_plot <- df %>%
          ggplot(aes(x = freq, y = reorder(`Gene symbol`, freq))) + 
            geom_col(fill="navy")
      }
    }
    
    basic_plot <- basic_plot +
      theme_bw() +
      theme(axis.text.y = element_text(size = 10)) +
      labs(
        y = "",
        x = "Gene frequency",
        title = paste0(
          "Genes with frequency in [",
          input$core_threshold[1], " - ", input$core_threshold[2], "] in ",
          input$selected_species
        )
      )
    
    basic_plot
  })
  
  return(ui)
}
