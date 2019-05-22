
###############################################################################
# update model comparison results with add variance explained for each model  #            #
# by fixed effects and phylogeny                                              #
# by traits and phylogeny                                                     #
# by traits alone                                                             #
#                                                                             #
# output the updated model comparison objects:                                #
# geographic_extent_model_comparison_disease.rData                            #
# host_range_model_comparisom_disease.rData                                   #
###############################################################################

rm(list=ls())

library(brms)
library(abind)

# load the function to calculate r2glmm
source("get_Nakagawa_R2glmm.R")

#for each of the impact responses 
for (i in 1:2){
  # load the model comparison object to bind the R2 results to
  load(paste0(wd, grep("model_comparison", list.files(), value = TRUE)[i])) 
  # load the null model required by Nakagawa methods for variance partitioning
  load(paste0("null_model_",
              paste0(strsplit(grep("model_comparison", 
                                   list.files(), 
                                   value = TRUE)[i], 
                              split = "_")[[1]][1:2], 
                     collapse = "_"),
              ".rData"))
  
  # get the model objects

  model_list <- grep("\\.csv", grep(paste0(strsplit(grep("model_comparison",
                                          list.files(),
                                          value = TRUE)[i], 
                                     split = "_")[[1]][1:2],
                            collapse = "_"),
                     list.files(list.dirs()[2]),
                     value = TRUE), invert = TRUE, value = TRUE)
  # reorder to match the order of the model index in the model comparison table
  model_list <- 
    paste0(list.dirs()[2],
           "/",
           gsub("[[:digit:]]|\\.rData", 
                "",  
                model_list), 
           1:length(model_list),
           ".rData")
                           
  # collect the data in an array for binding to the model results
  R2_results <- 
    array(dim = c(length(model_list),6,3),
          dimnames = list(NULL, 
                          c("R2_fixed_m",
                            "R2_traits_m",
                            "R2_fixed_c",
                            "R2_traits_c",
                            "ICC",
                            "ICCadj"),
                          c("est", "lower", "upper")))
  
  for(j in 1:length(model_list)){
    R2_results[j,,] <- 
      load_and_R2(model_name = model_list[j],
                  model_null = get(grep("null_model", ls(), value = TRUE)))
    print(j)
  }
  # bind it to the model results table
  assign(grep("model_comparison", 
              ls(), 
              value = TRUE),
         abind(get(grep("model_comparison", 
                        ls(), 
                        value = TRUE)), 
               R2_results, 
               along =2)
         )
  
  save(list =grep("model_comparison", 
                ls(), 
                value = TRUE),
       file = paste0(grep("model_comparison", 
                          ls(), 
                          value = TRUE), ".rData")
  )
  rm(list = c(grep("model_comparison", 
                 ls(), 
                 value = TRUE),
              grep("null_model", ls(), value = TRUE),
              "model_list"))
}

