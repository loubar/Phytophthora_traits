###############################################################################
# Functions to calculate variance components for fitted brms.fit models       #
# attributable to fixed and random effects (r2glmm, intra-class correlation)  #                                  
# following the procedure recommended by                                      #
# Nakagawa, S., Johnson, P., & Schielzeth, H. (2017).                         #
# The coefficient of determination R2 and intra-class correlation coefficient #
# from generalized linear mixed-effects models revisited and expanded.        #
# Journal of the Royal Society, Interface, 14(134)                            #
###############################################################################                                                                            #

# This function is based on the worked examples here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5636267/bin/rsif20170213supp2.pdf



get_R2_ICC <- function(model_f,        # fitted model with fixed effects (traits and years known)
                       model_null){    # null model with random effects only
  # Calculation of the variance in fitted values
  # Nakagawa example seems to be variance on the linear predictor scale
  # this value includes the variance explained by years known to science
  require(brms)
  
  sigma2_fixed <- data.frame(apply(fitted(model_f, 
                                          scale = "linear", # linear predictor scale
                                          robust = TRUE, 
                                          re_formula = NA, 
                                          summary = FALSE), 1, var))
  # it would be useful to pull out the variance explained by traits only
  # the code below tries to do this by fixing years known as a constant (the mean)
  # before calculating the variance: are we happy with this approach?
  sigma2_traits <- data.frame(apply(fitted(model_f, 
                                           scale = "linear", 
                                           robust = TRUE, 
                                           re_formula = NA, 
                                           summary = FALSE,
                                           newdata = data.frame(years_known = mean(model_f$data$years_known),
                                                                model_f$data[,grep("years_known",
                                                                                   colnames(model_f$data),
                                                                                   invert = TRUE)])), 1, var))
  # variance of phylogenetically structured random intercept
  sigma2_phylo <- posterior_samples(model_f, pars = "sd_species_name")^2
  sigma2_phylo_null <- posterior_samples(model_null, pars = "sd_species_name")^2
  
  # Intercept of the null model: expected grand mean on the linear predictor scale
  # accounting for phylogenetically structured random intercept
  # and observation-level random intercepts (additive dispersion)
  b0 <- as.numeric(fixef(model_null, summary = FALSE))
  
  # if the model is binomial (geographic extent) or bernoulli (disease symptoms)
  if(model_f$family$family == "binomial" ){
    
    # total variance on the linear predictor scale from the null model (see section 5 of paper: How to estimate lambda from data) 
    # sum of the phylogenetic variance (sd_species_name^2) and the additive dispersion parameter (sd_obs^2)
    sigma2_tau <- as.numeric((posterior_samples(model_null, pars = "^sd_species_name")^2) + (posterior_samples(model_null, pars = "^sd_obs")^2))
    
    # variance of observation-level random intercept to model overdispersion
    # (not used in the worked example) 
    # additive dispersion parameter
    sigma2_e <- posterior_samples(model_f, pars = "sd_obs")^2
    sigma2_e_null <- posterior_samples(model_null, pars = "sd_obs")^2
    
    n <- 159
    
    # convert to the response scale using Eqn. 6.8 (Nakagawa et al. 2017)
    pmean <- as.numeric(plogis(b0 - 0.5 * sigma2_tau * tanh(b0 * (1 + 2 * exp(-0.5 * sigma2_tau))/6)))
    
    # residual variance 
    # of null model
    sigma2_epsilon <- as.numeric(1/(n*pmean * (1 - pmean)))
  }
  if(model_f$family$family == "poisson"){ # see model 5 in the paper
    
    # total variance on the linear predictor scale from the null model (see section 5 of paper: How to estimate lambda from data) 
    # sum of the phylogenetic variance (sd_species_name^2) and the additive dispersion parameter (sd_obs^2)
    sigma2_tau <- as.numeric(posterior_samples(model_null, pars = "^sd_species_name")^2 + (posterior_samples(model_null, pars = "^sd_obs")^2))
    
    # variance of observation-level random intercept to model overdispersion
    # (not used in the worked example) 
    # additive dispersion parameter
    sigma2_e <- posterior_samples(model_f, pars = "sd_obs")^2
    sigma2_e_null <- posterior_samples(model_null, pars = "sd_obs")^2
    
    # estimate lambda: grand mean (= variance) via Eqn. 5.8 
    lambda <- exp(b0 + 0.5*sigma2_tau)
    # sigma2_epsilon is ln(1 + 1/lambda)
    sigma2_epsilon <- log(1 + 1/lambda)
  }
  if(model_f$family$family == "bernoulli" ){ 
    
    # total variance on the linear predictor scale from the null model (see section 5 of paper: How to estimate lambda from data) 
    # sum of the phylogenetic variance (sd_species_name^2) and the additive dispersion parameter (sd_obs^2)
    sigma2_tau <- as.numeric(posterior_samples(model_null, pars = "^sd_species_name")^2)
    
    # an additive overdispersion parameter doesn't make sense for binary models, so this parameter is null
    # no observation-level random effect in these models
    sigma2_e <- 0
    sigma2_e_null <- 0
    
    # convert to the response scale using Eqn. 6.8 (Nakagawa et al. 2017)
    pmean <- as.numeric(plogis(b0 - 0.5 * sigma2_tau * tanh(b0 * (1 + 2 * exp(-0.5 * sigma2_tau))/6)))
    
    # residual variance 
    # of null model
    sigma2_epsilon <- as.numeric(1/(pmean * (1 - pmean)))
  }
  #empty array for results (format for easy binding to model results table)
  results <-
    array(dim = c(1,6,3),
          dimnames = list(NULL, 
                          c("R2_fixed_m",
                            "R2_traits_m",
                            "R2_fixed_c",
                            "R2_traits_c",
                            "ICC",
                            "ICCadj"),
                          c("est", "lower", "upper")))
    
  # marginal (m) and conditional (c) R2
  # raw (ICC_phylo) and adjusted (ICCadj_phylo) intra-class correlation
  # return the summaries with the median estimate and lower and upper 95 credible intervals
  results[,"R2_fixed_m",]  <-  brms:::get_summary(robust = TRUE, sigma2_fixed/(sigma2_fixed + sigma2_phylo + sigma2_e + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_traits_m",] <-  brms:::get_summary(robust = TRUE, sigma2_traits / (sigma2_fixed + sigma2_phylo + sigma2_e + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_fixed_c",]  <-  brms:::get_summary(robust = TRUE, (sigma2_fixed + sigma2_phylo)/(sigma2_fixed + sigma2_phylo + sigma2_e + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_traits_c",] <-  brms:::get_summary(robust = TRUE, (sigma2_traits + sigma2_phylo)/(sigma2_fixed + sigma2_phylo + sigma2_e + sigma2_epsilon))[,c(1,3,4)]
  results[,"ICC",]         <-  brms:::get_summary(robust = TRUE, sigma2_phylo_null / (sigma2_phylo_null + sigma2_e_null + sigma2_epsilon))[,c(1,3,4)]
  results[,"ICCadj",]      <-  brms:::get_summary(robust = TRUE, sigma2_phylo / (sigma2_phylo + sigma2_e + sigma2_epsilon))[,c(1,3,4)]
  # note brms:::get_summary assumes version brms version 1.8.0
  # for more recent versions use brms::posterior_samples
  return(results)
  
}

load_and_R2 <- function(model_name, model_null){
  load(model_name)
  get_R2_ICC(model_f = model_fit,
             model_null)
}


