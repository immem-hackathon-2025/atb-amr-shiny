# plot candidate core genes for a selected species and core threshold
# take inputs from core_gene_species_controls.R

plotCoreGenesForSpecies <- function(core_threshold, selected_species) {
  # read in raw data, this has a merge of AFP and the species call + HQ
  ### TODO: REPLACE WITH INTERNAL DATA OBJECT
  ### IDEA: could remove filter to AMR, so we can plot the same info for virulence genes etc reported by AMRfp
  afp <-read_tsv("ATB_Enterobacter_AFP.tsv.gz") %>% filter(HQ) %>% filter(`Element subtype`=="AMR")
  
  # for a single species, plot candidate core genes
  # TODO: user settable
 
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
    ggplot(aes(x=freq, y=`Gene symbol`)) +
    geom_col()

}