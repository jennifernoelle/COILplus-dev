#' Function that updates the probability of occurence if not block-updated. 
#' Previously occurrence probabilities were specified via the arguments occur_B, occur_P
#' Now prior probabilities will be optionally supplied via the arguments p_occur_B, p_occur_P
#' We also add a new argument to the sampling vector to indicate whether or not these should be sampled
#' If not: we assume occurrence probabilities are known
#' If so: we assume occurrence probablities are known with uncertainty
#' 
#' 
#' Each species might occur in the study area or not. Even though probabilities
#' of occurence can be provided, the detected interactions also inform us of
#' whether the species occur. QUESTION: does this also hold for occurrence 
#' probabilities?
#' 
#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param occur Matrix. Rows are species, columns are studies. It includes the
#' prior probabilities of occurrence.
#' @param occor_others. Matrix Rows are species, columns are studies. It includes
#' the current values of occurrence for the other set of species, the ones that
#' are not being updated. 
#' @param probobs_curr Vector of current values for the probability of
#' observation for the set of species we are updating.
#' @param probobs_others Vector of current values for the probability of
#' observation for the other set of species, the ones that are not being
#' updated.
#' @param curr_inter Matrix with 0,1 entries for the current posterior samples
#' of possible interactions. The rows correspond to the species we are updating
#' and the columns correspond to the other set of species.
#' @param focus Array of three dimensions corresponding to the two sets of
#' species, and the different studies. Values are 0 or 1. This array represents
#' whether the current study would be willing to record the interactions of a
#' given species. The value should be 1 except for studies that are animal or
#' plant-oriented, where the value for species that are not of focus will be 0.


UpdOccurP <- function(occur_indicators, occur_prior_probs, curr_occur_probs, detected, mh_pprior_sd, mh_p_step) {
  
  num_obs <- nrow(occur_indicators)
  num_studies <- ncol(occur_indicators)
  
  # Return something like this to keep track of acceptance?
  new_value_p = curr_occur_probs
  accepted = matrix(0, nrow = num_obs, ncol = num_studies)
  
  # If the species were observed to interact in the study, the occurrence probability
  # is always equal to 1. First, I will draw the occurrence probabiilty
  # as if the species were not observed to interact. Then, I will assign the
  # ones observed interacting to 1.
  
  # Updating the occurrence probability as if the species was not detected in
  # the study. For all species, one study at a time.
  
  
  for (st in 1 : num_studies) {
    
    # Generating the proposal
    occur <- occur_indicators[, st]
    p_curr <- curr_occur_probs[, st]
    p_prior <- occur_prior_probs[, st]
    p_prop <- rtruncnorm(n = 1,  mean = p_curr, sd = mh_p_step, a = 0, b = 1)

    # Calculating the likelihood, prior, and proposal for proposed value for each proposed value: 
    log_prior_prop <- log(dtruncnorm(x = p_prop, a = 0, b = 1, mean = p_prior, sd = mh_pprior_sd))
    log_lik_prop <- log(p_prop^occur * (1-p_prop)^(1-occur))
    log_prop <- log(dtruncnorm(x = p_prop, a = 0, b= 1, mean = p_curr, sd = mh_pprior_sd))  
    
    # Calculating the likelihood, prior, and proposal for proposed value for each current value: 
    log_prior_curr <- log(dtruncnorm(x = p_curr, a = 0, b = 1, mean = p_prior, sd = mh_pprior_sd))
    log_lik_curr <- log(p_curr^occur * (1-p_curr)^(1-occur))
    log_curr <- log(dtruncnorm(x = p_curr, a = 0, b = 1, mean = p_prop, sd = mh_pprior_sd))
    
    AP <- log_lik_prop + log_prior_prop + log_curr
    AP <- AP - (log_lik_curr + log_prior_curr + log_prop)
    
    update <- log(runif(num_obs)) < AP
    new_value_p[update, st] <- p_prop[update]
    accepted[update, st] <- 1

  }
  
  # Now, if a species was detected to interact in a study, set their occurrence
  # probability to 1.
  wh_detected <- which(detected == 1)
  new_value_p[wh_detected] <- 1
  
  r <- list(new_value_p = new_value_p, accepted = accepted)
  return(r)
  
}

