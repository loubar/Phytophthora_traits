
###############################################################################
# to get the variance components, we need a null model with no fixed effects. #
#     	                                                                      #
# Nakagawa, S., Johnson, P., & Schielzeth, H. (2017).                         #
# The coefficient of determination R2 and intra-class correlation coefficient #
# from generalized linear mixed-effects models revisited and expanded.        #
# Journal of the Royal Society, Interface, 14(134)                            #
#                                                                             #
# our null models in the model comparison contain years known to science      #
# so we'll need to fit null models now with only the random effects           #
#                                                                             #
# outputs:                                                                    #
# null_model_geographic_extent.rData                                          #
# null_model_host_range.rData                                                 #
###############################################################################

rm(list=ls())

library(brms)

load("Phytophthora_covariance.rData")
load("disease_traits/impact_geographic_extent_disease_1.rData")
# load an example fitted model object 
# use the trait data from these saved objects to ensure the null models 
# are fitted to the same subset of species

# set priors
prior <- c(prior(normal(0, 50), "Intercept"),
           prior(student_t(4, 0, 1), "sd")) 


# fit the null models with only random effects
# our null models include fixed effect of years_known_to_science


# binomial model with additive dispersion
null_model_geographic_extent <-
  brm(impact_geographic_extent | trials(159) ~ 1 + (1 | species_name) + (1 | obs), 
      data = model_fit$data, 
      family = binomial(link = "logit"),
      cov_ranef = list(species_name = Phytophthora_covariance),
      prior = prior,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 15),
      cores = 4,
      iter = 5000,
      warmup = 4000,
      save_dso = TRUE)
save(null_model_geographic_extent,
     file = "null_model_geographic_extent.rData")




# poisson model with additive dispersion
null_model_host_range <-
  brm(impact_host_range ~ 1 + (1 | species_name) + (1 | obs), 
      data = model_fit$data, 
      family = poisson(link = "log"),
      cov_ranef = list(species_name = Phytophthora_covariance),
      prior = prior,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 15),
      cores = 4,
      iter = 5000,
      warmup = 4000,
      save_dso = TRUE)
save(null_model_host_range,
     file = "null_model_host_range.rData")

