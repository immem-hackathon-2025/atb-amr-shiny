# Visualization of AMRfinderplus results from AllTheBacteria

A Shiny app developped at the 12th Microbial Bioinformatics Hackathon (IMMEM XIV)!


## AllTheBacteria Preprocessing

AllTheBacteria's AMRfinderplus results were preprocessed to make a smaller sized input file that is more manageable for the Shiny app.

AllTheBacteria's AMRfinderplus results and species calls were merged and then filtered by the following columns to keep:
- HQ = TRUE (keep only high quality species calls)
- status = PASS (keep genomes where AMRFP sucessfully completed)

The remaining columns that are kept as a simplified/preprocessed input file:
- Name, Gene symbol, Hierarchy node, Element type, Element subtype, % Coverage of reference sequence, % Identity to reference sequence, Method, Class, Subclass

ATB_preprocessing_perspecies.Rmd - R code to preprocess AMRfp results for a single, specified species

ATB_preprocessing.Rmd - R code to preprocess AMRfp results for all species
