
############################################################################################
# Process the raw trait data and phylogenetic tree
# fit pure latent variable models and a null model for comparison
# compile a list of candidate models for multi-model model selection
# output the following objects
#  a) latent variable model         : fit_LV_model_disease.rData
#  b) null LV model                 : fit_null_model_disease.rData
#  c) the latent variable estimates : lvs_disease.rData
#  d) phylogenetic covariance matrix: Phytophthora_covariance.rData
#  e) processed trait dataframe     : traits.rData
#  f) list of candidate models      : all_predictors_disease.rData 
############################################################################################
 


rm(list=ls())

library(R2jags)
library(ggplot2)
library(cluster)
library(MuMIn)
library(ape)
library(RNeXML)
library(brms)
library(boral)
library(phytools)

set.seed(123)



############################################################################################
# read in the host range, disease and trait data
############################################################################################


# host range data
host_range <- read.csv("data/2018-04-18_impact_host_range.csv",
                          stringsAsFactors = FALSE)

# disease symptom data
disease_symptoms <- read.csv("data/2019-01-28_disease_symptoms.csv",
                            stringsAsFactors = FALSE)[,1:3]

# trait database
traits <- read.csv("data/2018-01-29_Phytophthora_trait_database_merged.csv",
                   stringsAsFactors = FALSE,
                   na.strings = c("NA", "#N/A"))

############################################################################################
# read in and process the geographic extent data
############################################################################################

extent_db <- read.csv("data/2018-04-17_country_level_occurrence_database.csv", 
                      stringsAsFactors = FALSE)

extent_db <- extent_db[!is.na(extent_db$spname),]

# Phytophthora alni records are those where the subspecies 
# was not identified: remove these
# where the subspecies of hybrid 
# 
extent_db <- extent_db[!extent_db$spname %in% "Phytophthora alni" ,]

# now remeove the records that refer to eradicated, extinct or absent
#"Alert List",
#"Absen",
#"Phytosanitary incident",
#"Extinct"
# Transient
# "Eradicated"
extent_db <- 
  extent_db[
    grep("Alert List|Absen|Phytosanitary incident|Extinct|eradicat|transient",
         extent_db$Introduction_code,
         invert = TRUE),]



species_by_country <- 
  acast(data = extent_db, 
        spname~country_clean_std, 
        length) 
# convert number of records to presence absence
species_by_country[species_by_country > 0] <- 1

geographic_extent <- data.frame(species_name = names(rowSums(species_by_country)), 
                                impact_geographic_extent = rowSums(species_by_country),
                                row.names = NULL)

# how many species do we have geographic extent data for?
nrow(geographic_extent)

# clean the species names to match the trait database
setdiff(geographic_extent$species_name, traits$species_name)
geographic_extent$species_name  <- gsub("?", "x ", geographic_extent$species_name)

# merge each of these with the trait database
# check species names match first
all(geographic_extent$species_name %in% traits$species_name)
all(host_range$species_name %in% traits$species_name)
all(disease_symptoms$species_name %in% traits$species_name)

traits <- 
  merge(traits, merge(geographic_extent, merge(disease_symptoms, host_range)))
# 1 + years known to science 
traits$years_known <- as.integer(substr(Sys.time(), 1, 4)) - traits$date
# number of countries sampled for binomial model of geographic extent
traits$n_trials <- length(unique(extent_db$country_clean_std))


# choose the traits to model impact with
included_traits <- c("proliferation",
                     "caduceus",
                     "oospore_wall_index",
                     "chlamydospores",
                     "hyphal_swelling",
                     "root_disease",
                     "foliar_disease",
                     "oospores",
                     "gr_at_opt",
                     "temp_opt")

impact_metrics <- c("impact_geographic_extent",
                    "impact_host_range")



cont_variables <- 
  names(which(sapply(traits[,included_traits], is.double)))

binary_variables <- 
  c(names(which(sapply(traits[,included_traits], is.character))), "root_disease", "foliar_disease")

# set the oospore wall index for species without oospores (sterile) to the mean value
traits$oospore_wall_index[which(traits$ho_he_s == "S")] <- mean(traits$oospore_wall_index, 
                                                                na.rm = TRUE) 



# convert traits to binary variables for binomial modelling
for(i in grep("disease",
                       binary_variables,
                       value = TRUE,
                       invert=TRUE)){
  traits[,i] <- ifelse(grepl("^N", traits[,i]), 
                                       yes = 0,
                                       no = 1)
}

# the apprropriate subset of data is extracted in the model fitting functions on the cluster

# how many species have complete trait data for host range models (including disease traits?)
nrow(traits[complete.cases(traits[,c(included_traits, "years_known", "impact_host_range")]),])
#without disease symptoms
nrow(traits[complete.cases(traits[,c(grep("disease",included_traits, invert = TRUE, value = TRUE), "years_known", "impact_host_range")]),])

# how many species have complete trait data for geographic extent models
nrow(traits[complete.cases(traits[,c(included_traits, "years_known", "impact_geographic_extent")]),])
#without disease symptoms
nrow(traits[complete.cases(traits[,c(grep("disease",included_traits, invert = TRUE, value = TRUE), "years_known", "impact_geographic_extent")]),])

#nrow(traits_no_disease <- 
#       traits[complete.cases(traits[,grep("disease", 
#                                          included_traits, 
#                                          value = TRUE, 
#                                          invert = TRUE)]),])

# we could use the predict function in brms to sample from the posterior distribution  to
# and to predict their host ranges and geographic extents
# into the future using the parameters from our fitted models,
# their position in the phylogeny and their trait values, conditional on years known
# to science.  This would be a nice excercise to finish the paper with showing how the early-
# warning systemn could be applied predictively
# it would be fun to do this for recently described species and/or those without impact data

# which species without impact data could we predict for?
x <- traits[complete.cases(traits[, included_traits]),]
x[is.na(x$impact_geographic_extent),"species_name"]
traits[is.na(traits$impact_geographic_extent),"species_name"]




##############################################################################################
# fit the latent variable models and add the latent variable estimates to the trait database # 
##############################################################################################
nLVs <- 2

fit_LV_model_disease <- 
  boral(y = traits[,included_traits],
        lv.control = list(num.lv = nLVs),
        family = ifelse(colnames(traits[,included_traits]) %in% binary_variables, 
                        yes = "binomial", 
                        no = "normal"),
        model.name = "data/LV_model_disease.txt",
        save.model = TRUE
        )

fit_null_model_disease <- 
  boral(y = traits[,included_traits],
        lv.control = list(num.lv = 0),
        family = ifelse(colnames(traits[,included_traits]) %in% binary_variables, 
                        yes = "binomial", 
                        no = "normal"),
        model.name = "data/null_model_disease.txt",
        save.model = TRUE
  )
# deviance explained
(fit_null_model_disease$jags.model$BUGSoutput$median$deviance-fit_LV_model_disease$jags.model$BUGSoutput$median$deviance)/fit_null_model_disease$jags.model$BUGSoutput$median$deviance


plot(fit_LV_model_disease)

lvs_disease <- as.data.frame(fit_LV_model_disease$lv.median[,1:nLVs])
#lvs_no_disease <- as.data.frame(fit_LV_model_no_disease$lv.median[,1:nLVs])



save(fit_LV_model_disease,
     file = "data/fit_LV_model_disease.rData")
save(fit_null_model_disease,
     file = "data/fit_null_model_disease.rData")
save(lvs_disease,
     file = "data/lvs_disease.rData")




############################################################################################
# prepare the phylogenetic covariance matrix
############################################################################################
Phytoph_tree <-
  read.tree("data/updated - Posterior output.newick")

# clean the tip labels to match the trait database
Phytoph_tree$tip.label <- gsub("Phytophthora_", "Phytophthora ", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("x_", "x ", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("austrocedrae",
                               "austrocedri",
                               Phytoph_tree$tip.label)
Phytoph_tree$tip.label[grep(" sp", Phytoph_tree$tip.label)] <- 
  sapply(grep(" sp", Phytoph_tree$tip.label),
         function(i)
           paste(strsplit(Phytoph_tree$tip.label, split = " sp")[[i]][1], 
                 paste0("'", strsplit(Phytoph_tree$tip.label, split = " sp")[[i]][2], "'"), 
                 collapse = " "))
Phytoph_tree$tip.label <- gsub("gondwanense", "gondwanensis", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("_", "", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("emzansi", "emanzi", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("bisheria", "bishii", Phytoph_tree$tip.label)

# check all species in the trait database are represented in the phylogeny
setdiff(traits$species_name, Phytoph_tree$tip.label)

Phytoph_tree <- multi2di(Phytoph_tree)
Phytoph_tree$edge.length[Phytoph_tree$edge.length==0] <- max(nodeHeights(Phytoph_tree))*0.001

Phytoph_tree <-
  chronopl(phy = Phytoph_tree,
           lambda = 1)
is.ultrametric(Phytoph_tree)

inv.phylo <- 
  MCMCglmm::inverseA(Phytoph_tree, nodes = "TIPS", scale = TRUE)
Phytophthora_covariance <- solve(inv.phylo$Ainv)
rownames(Phytophthora_covariance) <- rownames(inv.phylo$Ainv)

save(Phytophthora_covariance,
     file = "data/Phytophthora_covariance.rData")


#################################################################################
# prepare the Martin et al. multi-gene phylogeny also for comparison of results #
#################################################################################

# Phytoph_phylogeny <-
#   nexml_read(x = "http://purl.org/phylo/treebase/phylows/tree/TB2:Tr65776?format=nexml")   
# 
# # coerce the NeXML object to a phylo object recognised by
# # ape, ade4, phylobase, adephylo
# Phytoph_tree <- get_trees(Phytoph_phylogeny)
# 
# # prune away one of the citrophthora isolates as duplicated tip labels won't work when matching to the trait data
# Phytoph_tree <- drop.tip(phy = Phytoph_tree, tip = c("Phytophthora citrophthora P10341", "Phytophthora palmivora arecae P10213"))
# 
# # remove isolate numbers
# Phytoph_tree$tip.label <- gsub(" P$", "", gsub("[[:digit:]]", "", Phytoph_tree$tip.label))
# 
# 
# 
# # clean the tip labels to match the trait database
# #Phytoph_tree$tip.label <- gsub("Phytophthora ", "P. ", Phytoph_tree$tip.label)
# 
# Phytoph_tree$tip.label <- gsub("alni", "x alni", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("andina", "x andina", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("katsurae", "castaneae", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("kelmania", "'kelmania'", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("lagoariana", "'lagoariana'", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("ohioensis", "'ohioensis'", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("personii", "'personii'", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("palmivora arecae", "palmivora", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("cinnamomi parvispora", "parvispora", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("salixsoil", "lacustris", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("sansomea", "sansomeana", Phytoph_tree$tip.label)
# Phytoph_tree$tip.label <- gsub("bisheria", "bishii", Phytoph_tree$tip.label)
# 
# is.rooted(Phytoph_tree) # FALSE
# #There are two outgroups, *Pythium undulata* and *Pythium vexans*, 
# # forming a polytomy at the node with the Phytophthora genus.  
# # If we drop the *Pythium vexans* tip, the tree should be properly rooted.
# Phytoph_tree <- drop.tip(phy = Phytoph_tree, tip = "Pythium vexans")
# 
# 
# is.ultrametric(Phytoph_tree) # FALSE
# # No, we'll have to estimate the ages of the nodes by assuming branch lengths 
# # are proportional to the number of substitutions.  Is this a reasonable approach?  
# # Most of the trait evolution tools assume the phylogeny is ultrametric.
# Phytoph_tree <- chronopl(phy = Phytoph_tree,
#                          lambda = 1)
# 
# 
# is.rooted(Phytoph_tree) # TRUE
# is.ultrametric(Phytoph_tree) # TRUE
# 
# setdiff(traits$species_name, Phytoph_tree$tip.label)
# 
# # prepare the covariance matrix
# inv.phylo <- 
#   MCMCglmm::inverseA(Phytoph_tree, nodes = "TIPS", scale = TRUE)
# Phytophthora_covariance <- solve(inv.phylo$Ainv)
# rownames(Phytophthora_covariance) <- rownames(inv.phylo$Ainv)
# 
# # how many species can we include in the models using multi-gene phylogeny?
# length(intersect(Phytoph_tree$tip.label, traits$species_name))
# # 78
# 
# save(Phytophthora_covariance,
#      file = "Phytophthora_covariance_multigene.rData")

############################################################################################
# set up a list of candidate models for model comparison                                   #
############################################################################################

# For the individual trait models
# rescale continuous variables by 2 SDs as per Gelman 
# so they have variance 
# approximately the same as the binary variables 
# this means the parameter estimates will be comparable between binary and continuous variables
# we don't centre the trait variables as this would mean the binary variables 
# would become -0.5 and + 0.5.  Boral needs binary varables to be coded 0 and 1
# the parameter estimates in the brms models can then be interpreted as the number of units increase
# in impact (hosts, countries) with a 2SD increase in trait x

scale2sd <- function(var){
  return(var/(2*sd(var, na.rm = TRUE)))
}
traits[,c(cont_variables, "years_known")] <- 
  sapply(traits[,c(cont_variables, "years_known")], scale2sd)

# set up the global model predictors
predictors <- list()
predictors[[1]] <- c("years_known", included_traits, "foliar_disease:root_disease")
predictors[[2]] <- c("years_known", "lv1", "lv2", "lv1:lv2")
# use dredge to get a list of all possible candidate models (subsets of predictors) from the global models with each set of predictors (individual traits, clusters and axes of trait-space)



model_sets <- list()
dummy_traits <- as.data.frame(sapply(1:length(c("impact", included_traits, "years_known")),
                                     function(i) 
                                       rnorm(n = 20)))
colnames(dummy_traits) <- c("impact", included_traits, "years_known")
model_traits <- 
    lm(formula = as.formula(
      paste("impact",
            "~", 
            paste(predictors[[1]], 
                  collapse = " + "), 
            collapse = " ")),
      data = dummy_traits,
      na.action = "na.fail")
dredge_global_model <- 
    dredge(global.model = model_traits,
           evaluate = FALSE)
  
  
  model_sets[[1]] <-
    #    model_sets <-
    gsub(" \\+ 1", 
         "", 
         sapply(1:length(dredge_global_model),
                function(x) 
                  strsplit(split = "\\ ~ ",
                           as.character(dredge_global_model[[x]]))[[2]][2]))

  dummy_traits <- as.data.frame(sapply(1:length(c("years_known","impact", "lv1", "lv2")),
                                       function(i) 
                                         rnorm(n = 20)))
  colnames(dummy_traits) <- c("years_known","impact", "lv1", "lv2")

  model_LVs <- 
    lm(formula = as.formula(
      paste("impact",
            "~", 
            paste(predictors[[2]], 
                  collapse = " + "), 
            collapse = " ")),
      data = dummy_traits,
      na.action = "na.fail")
  dredge_global_model <- 
    dredge(global.model = model_LVs,
           evaluate = FALSE)
  
  
  model_sets[[2]] <-
    gsub(" \\+ 1", 
         "", 
         sapply(1:length(dredge_global_model),
                function(x) 
                  strsplit(split = "\\ ~ ",
                           as.character(dredge_global_model[[x]]))[[2]][2]))
  
  
# remove the intercept only model: we're going to fit these manually 
# so they can be used a benchmark to estimate variance explained by each model
all_predictors <-
  unique(unlist(model_sets))
#all_predictors <- all_predictors[grep("^1$", 
#                                      all_predictors,
#                                      invert = TRUE)]

# how many models to fit for model selection (with and without disease symptoms)
length(all_predictors_disease <- all_predictors)
# length(all_predictors_no_disease <- grep("disease", 
#                                        all_predictors, 
#                                        value = TRUE, 
#                                        invert = TRUE))
all_predictors_disease <- grep("years_known", 
                               all_predictors_disease, 
                               value = TRUE, 
                               invert = TRUE)
length(all_predictors_disease)


save(traits,
     file = "data/traits.rData")
save(all_predictors_disease,
     file = "data/all_predictors_disease.rData")



