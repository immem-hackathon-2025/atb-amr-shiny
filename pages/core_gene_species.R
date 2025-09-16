# plot gene frequencies for a selected species and frequency range
coreGenesForSpeciesPage <- function(afp, input, output) {

  
  # --- Build species list from distinct samples (Name) per Species ---
  # Only include species with at least 10 unique samples
  species_counts <- afp %>%
    distinct(Name, Species) %>%                  # unique sample-species pairs
    group_by(Species) %>%
    summarise(nspp = n(), .groups = "drop") %>%
    filter(nspp >= 10) %>%
    arrange(desc(nspp)) %>%
    collect()
  
  species_choices <- setNames(
    species_counts$Species,
    paste0(species_counts$Species, " (n=", species_counts$nspp, ")")
  )
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "selected_species",
          "Choose a species to explore its gene frequency:",
          choices = species_choices,
          selected = if (length(species_choices)) species_counts$Species[[1]] else NULL
        ),
        # select gene threshold for core_gene_species plot
        sliderInput(
          "core_threshold",
          "Select a frequency range for genes to include:",
          min = 0, max = 1, value = c(0.6, 1.0)
        ),
        # sliders are proportions (0..1); converted below to percent to match columns
        sliderInput(
          "identity_threshold",
          "Minimum nucleotide identity (proportion):",
          min = 0.5, max = 1.0, value = 0.9
        ),
        sliderInput(
          "coverage_threshold",
          "Minimum coverage (proportion):",
          min = 0.5, max = 1.0, value = 0.9
        )
      ),
      mainPanel(
        plotOutput("coreGeneSpeciesPlot", height = "calc(100vh - 200px)"),
        shiny::uiOutput("coreGeneSpeciesDownloadButton")
      )
    )
  )
  
  # Reactive data for selected species and thresholds
  coreGeneSpecies <- reactive({
    req(input$selected_species, input$core_threshold, input$identity_threshold, input$coverage_threshold)
    
    # Convert proportions [0,1] to percentages [0,100] to match the AFP columns
    id_min  <- input$identity_threshold  * 100
    cov_min <- input$coverage_threshold  * 100
    
    # Count how many distinct samples contain each gene (within species)
    # Then compute frequency = n / n_samples_for_species
    # Push everything to DuckDB, collect only the filtered result.
    # Column names with spaces or % need backticks.
    spp_tbl <- afp %>% filter(Species == !!input$selected_species)
    
    n_samples <- spp_tbl %>%
      summarise(n = n_distinct(Name)) %>%
      collect() %>%
      pull(n)
    
    validate(need(n_samples > 0, paste0("No rows for species: ", input$selected_species)))
    
    gene_freq_tbl <- spp_tbl %>%
      filter(`% Coverage of reference sequence` >= !!cov_min) %>%
      filter(`% Identity to reference sequence`  >= !!id_min) %>%
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
  output$coreGeneSpeciesDownloadButton <- renderUI({
    # IconButton is assumed to be defined elsewhere in your app
    IconButton("downloadCoreGeneSpecies", "data_dl", "Download")
  })
  
  output$downloadCoreGeneSpecies <- downloadHandler(
    filename = function() {
      paste0("core_gene_species_", gsub("\\s+", "_", input$selected_species), ".tsv")
    },
    content = function(file) {
      df <- coreGeneSpecies()
      write_tsv(df, file)
    }
  )
  
  # Plot
  output$coreGeneSpeciesPlot <- renderPlot({
    df <- coreGeneSpecies()
    
    validate(need(nrow(df) > 0, "No genes in the selected frequency range with current thresholds."))
    
    ggplot(df, aes(x = freq, y = reorder(`Gene symbol`, freq), fill = Class)) +
      geom_col() +
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
  })
  
  return(ui)
}
