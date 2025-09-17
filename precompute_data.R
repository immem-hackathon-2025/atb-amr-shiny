#!/usr/bin/env Rscript
# Precompute species, genus, and gene data for each element type
# Run this script once to generate cached data files

library(shiny)
library(tidyverse)
library(dplyr)
library(DBI)
library(duckdb)

message("Starting data precomputation...")

# Connect to database
con <- dbConnect(duckdb())

# Point to your Parquet dataset (Hive-style partitioned by Genus)
afp_amr <- tbl(con, "read_parquet('data/amr_by_genus/**/*.parquet', hive_partitioning = true)")
afp_stress <- tbl(con, "read_parquet('data/stress_by_genus/**/*.parquet', hive_partitioning = true)")
afp_virulence <- tbl(con, "read_parquet('data/virulence_by_genus/**/*.parquet', hive_partitioning = true)")

connections <- list(
  "AMR" = afp_amr,
  "STRESS" = afp_stress,
  "VIRULENCE" = afp_virulence
)

# Create cache directory if it doesn't exist
cache_dir <- "data"
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir)
  message("Created cache directory")
}

# Precompute species, genus, and gene data for each connection
message("Precomputing species, genus, and gene data for each element type...")

species_data <- list()
genus_data <- list()
gene_data <- list()

for (element_type in names(connections)) {
  message(paste("Computing species data for", element_type, "..."))
  
  # Calculate species counts (only include species with at least 10 unique samples)
  species_counts <- connections[[element_type]] %>%
    distinct(Name, Species) %>%                  # unique sample-species pairs
    group_by(Species) %>%
    summarise(nspp = n(), .groups = "drop") %>%
    filter(nspp >= 10) %>%
    arrange(desc(nspp)) %>%
    collect()
  
  # Create choices for selectize input
  species_choices <- setNames(
    species_counts$Species,
    paste0(species_counts$Species, " (n=", species_counts$nspp, ")")
  )
  
  # Store both counts and choices
  species_data[[element_type]] <- list(
    counts = species_counts,
    choices = species_choices
  )
  
  message(paste("Computing genus data for", element_type, "..."))
  
  # Add genus extraction and compute genus data
  afp_with_genus <- connections[[element_type]] %>%
    mutate(
      Genus = sql("substr(Species, 1, instr(Species, ' ') - 1)"),
      species_name = sql("substr(Species, instr(Species, ' ') + 1)")
    )
  
  # Calculate genus counts
  genus_counts_base <- afp_with_genus %>%
    distinct(Name, Genus) %>%                  # unique sample-genus pairs
    group_by(Genus) %>%
    summarise(n = n(), .groups = "drop") 
  
  genus_counts_final <- afp_with_genus %>%
    distinct(Species, Genus) %>%                  # unique species-genus pairs
    group_by(Genus) %>%
    summarise(nspp = n(), .groups = "drop") %>%
    filter(nspp >= 10) %>%
    left_join(genus_counts_base, by="Genus") %>%
    arrange(desc(nspp)) %>% 
    collect()
  
  # Create genus choices
  genus_choices <- genus_counts_final$Genus
  names(genus_choices) <- paste(
    genus_counts_final$Genus,
    " (n=", genus_counts_final$n,
    " in ", genus_counts_final$nspp,
    " species)",
    sep = ""
  )
  
  # Store genus data
  genus_data[[element_type]] <- list(
    counts = genus_counts_final,
    choices = genus_choices
  )
  
  message(paste("Computing gene data for", element_type, "..."))
  
  # Calculate gene counts
  gene_counts <- connections[[element_type]] %>%
    group_by(`Gene symbol`) %>%
    summarise(
      n = n(),
      nspp = n_distinct(Species)
    ) %>%
    arrange(desc(n)) %>%
    filter(!is.na(`Gene symbol`)) %>%
    collect()
  
  # Create gene choices
  gene_choices <- gene_counts$`Gene symbol`
  names(gene_choices) <- paste(
    gene_counts$`Gene symbol`,
    " (n=", gene_counts$n,
    " in ", gene_counts$nspp,
    " species)",
    sep = ""
  )
  
  # Store gene data
  gene_data[[element_type]] <- list(
    counts = gene_counts,
    choices = gene_choices
  )
}

# Save precomputed data to disk
message("Saving precomputed data to cache files...")

saveRDS(species_data, file.path(cache_dir, "species_data.rds"))
saveRDS(genus_data, file.path(cache_dir, "genus_data.rds"))
saveRDS(gene_data, file.path(cache_dir, "gene_data.rds"))

# Save metadata about when this was computed
metadata <- list(
  computed_at = Sys.time(),
  element_types = names(connections),
  species_counts = sapply(species_data, function(x) nrow(x$counts)),
  genus_counts = sapply(genus_data, function(x) nrow(x$counts)),
  gene_counts = sapply(gene_data, function(x) nrow(x$counts))
)

saveRDS(metadata, file.path(cache_dir, "metadata.rds"))

# Close database connection
dbDisconnect(con)

message("Data precomputation complete!")
message("Cache files saved:")
message(paste("  -", file.path(cache_dir, "species_data.rds")))
message(paste("  -", file.path(cache_dir, "genus_data.rds")))
message(paste("  -", file.path(cache_dir, "gene_data.rds")))
message(paste("  -", file.path(cache_dir, "metadata.rds")))

# Print summary
message("\nSummary:")
for (et in names(connections)) {
  message(paste0("  ", et, ":"))
  message(paste0("    Species: ", nrow(species_data[[et]]$counts)))
  message(paste0("    Genera: ", nrow(genus_data[[et]]$counts)))
  message(paste0("    Genes: ", nrow(gene_data[[et]]$counts)))
}
