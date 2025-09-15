plotCoreGenesAcrossGenus <- function(afp=afp, core_threshold, selected_genus, min_genomes_per_species) {
  
  afp <- afp %>% filter(grepl(selected_genus, Species))
  
  # total number per species
  n_per_species <- afp %>% 
    group_by(Name, Species) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(Species) %>% 
    summarise(nspp=n()) %>% 
    arrange(-nspp)
  
  ### TODO: allow user to select node instead of gene, as the unit of measurement
  # gene frequency per species
  gene_count_per_spp <- afp %>% 
    group_by(Name, `Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(`Gene symbol`, Class, Subclass, Species, `Element type`) %>% 
    count() %>% 
    left_join(n_per_species, by="Species") %>% 
    mutate(freq=n/nspp)
  
  # extract core genes and display freq per species
  gene_count_per_spp %>% filter(nspp>min_genomes_per_species & freq>core_threshold) %>%
    mutate(label=paste0(Species, " (n=", nspp, ")")) %>%
    arrange(-nspp) %>%
    ggplot(aes(y=label, x=freq)) + 
    geom_col(position='dodge') + 
    facet_wrap(~`Gene symbol`) + 
    theme(legend.position="bottom")
}