# Phytophthora_traits
## Trait-based approaches for plant health risk assessment for *Phytophthora* pathogens

The data and code reproduce the analyses in:

1) Barwell *et al.* (in prep.) A conceptual trait-based framework for invasion risk assessment of oomycete plant pathogens

2) Barwell *et al.* (in prep.), Trait-based approaches for predicting future global impacts in the genus *Phytophthora*

These studies develop and test ideas for a trait-based early warning system for pathogens for use in plant health risk assessment.

The data folder contains the *Phytophthora* trait database, ITS phylogeny, host range and country-level distribution data and habitat preferences. There is also a file of supporting information describing the data in the trait database.   

The file explore_trait_database.RMD reproduces the analyses in the main text of Barwell et al. (in prep.) A conceptual trait-based framework for invasion risk assessment of oomycete plant pathogensthe paper.  It generates the explore_trait_database.html document to help users visualise trait evolution and patterns of trait covariance.  

The numbered scripts reproduce the analyses in Barwell et al. (in prep.), Trait-based approaches for predicting future global impacts in the genus *Phytophthora* . The scripts are numbered in the order they should be run. Each script outputs the R objects required for the next stage of the analysis.

As the full set of candidate models includes 1284 models, the model comparison was run in parallel on a cluster. Reproducing step 2 (2_fit_global_impact_models_disease.R) could, therefore, take a very long time. To enable users to complete the workflow from start to finish in a reasonable amount of time, the outputs from 3_collate_model_results.R are also provided, meaning the model fitting (steps 2 and 3) can be skipped if required.
