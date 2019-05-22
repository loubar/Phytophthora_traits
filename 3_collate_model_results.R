
###########################################################################
# collate the parameter estimates and information criteria for            # 
# for all 1284 models ready for ranking                                   #
# output a compiled model results file for each impact response variable: #
#                                                                         #
# geographic_extent_model_comparison_disease.rData                        #
# host_range_model_comparison_disease.rData                               #
###########################################################################
rm(list = ls())
load("data/all_predictors_disease.rData")

# first check all the models have been fitted - if all models are present, both host_missing and geo_missing should be length 0
# If not, go back and rerun 3_fit_global_impact_models_disease.R
# It will automatically check which models have results already and only rereun the missing models
# Sometimes a node throws up an error I can't find the root of and the job is exited before all models have been fitted
# I think snowfall is a bit temperamental with any functions that have a built-in "cores = " argument, even if it isn't used
length(host_missing <- setdiff(1:length(all_predictors_disease), 
			as.numeric(gsub("\\D", 
				   "", 
				   grep("host_range", 
					grep("\\.rData", 
					     list.files("disease_traits"), 
					     value = TRUE), 
				   	value = TRUE)))))

length(geo_missing <- setdiff(1:length(all_predictors_disease), 
			as.numeric(gsub("\\D", 
				   "", 
				   grep("geographic_extent", 
					grep("\\.rData", 
					     list.files("disease_traits"), 
					     value = TRUE), 
				   	value = TRUE)))))
 


#create an empty array for host range and geographic extent to populate with the point estimates and lower and upper 95% credible intervals of the 
# parameter estimates
geographic_extent_model_comparison_disease <-
	array(dim = c(length(all_predictors_disease),
		      length(grep("object_name", colnames(read.csv("disease_traits/est_impact_geographic_extent_disease_1.csv")), value = TRUE, invert = TRUE)
),
		      3),
               dimnames = list(rep(NULL,length(all_predictors_disease)) , grep("object_name", colnames(read.csv("disease_traits/est_impact_geographic_extent_disease_1.csv")), value = TRUE, invert = TRUE)
, c("est", "lower", "upper")))
host_range_model_comparison_disease <- geographic_extent_model_comparison_disease

########################### geographic extent  #######################


for (i in c("est", "lower", "upper")){
	x <- do.call(rbind, lapply(paste0("disease_traits/", grep(paste0(i, "_impact_geographic_extent_"), 
                      grep("\\.csv", 
                           list.files("disease_traits"), 
                           value = TRUE), 
                      value = TRUE)), 
                 read.csv,
		 stringsAsFactors = FALSE))
	z <- x[order(x$model_index),]
	z$object_name <- NULL
	geographic_extent_model_comparison_disease[,,i] <- as.matrix(z)[,dimnames(geographic_extent_model_comparison_disease)[[2]]]
}

save(geographic_extent_model_comparison_disease,
     file = "geographic_extent_model_comparison_disease.rData")


############################ host range ###############################

for (i in c("est", "lower", "upper")){
	x <- do.call(rbind, lapply(paste0("disease_traits/", grep(paste0(i, "_impact_host_range_"), 
                      grep("\\.csv", 
                           list.files("disease_traits"), 
                           value = TRUE), 
                      value = TRUE)), 
                 read.csv,
		 stringsAsFactors = FALSE))
	z <- x[order(x$model_index),]
	z$object_name <- NULL
	host_range_model_comparison_disease[,,i] <- as.matrix(z)[,dimnames(host_range_model_comparison_disease)[[2]]]
}

save(host_range_model_comparison_disease,
     file = "host_range_model_comparison_disease.rData")



