# plot candidate core genes for a selected species and core threshold
# take inputs from core_gene_species_controls.R

plotCoreGenesForSpecies <- function(afp, core_threshold, selected_species) {
 
  afp_this_spp <- afp %>% filter(Species==selected_species)
  n_this_spp <- length(unique(afp_this_spp$Name))
  
  # TODO: make pretty, annotate with class/subclass
  afp_this_spp %>%
    group_by(Name, `Gene symbol`, Class, Subclass) %>% 
    count() %>% distinct() %>% ungroup() %>% 
    group_by(`Gene symbol`, Class, Subclass) %>% 
    count() %>% 
    mutate(freq=n/length(unique(afp_this_spp$Name))) %>%
    filter(freq>=core_threshold) %>% 
    ggplot(aes(x=freq, y=`Gene symbol`, fill=Class)) +
    geom_col() + 
    theme(axis.text.y=element_text(size=10)) + 
    labs(y="", x="Gene frequency", title=paste0("Genes with freq >= ",core_threshold," in ", selected_species))

}