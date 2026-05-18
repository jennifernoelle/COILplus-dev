# Blocked sampler using
# Occurrence and probability are sampled independently and then jointly accepted or rejected


#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param occur_prior_probs Matrix.  Rows are species, columns are studies. 
#' It includes the user provided prior probabilities of occurrence
#' @param p_curr (Previously 'occur'). Matrix. Rows are species, columns are studies. 
#' It includes the probabilities of occurrence at the current step
#' @param occur_others Matrix Rows are species, columns are studies. It includes
#' the current values of occurrence for the other set of species, the ones that
#' are not being updated. 
#' @param probobs_curr Vector of current values for the probability of detection/
#' observation for the set of species we are updating.
#' @param probobs_others Vector of current values for the probability of detection/
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
#' 
#' 

UpdOccurP_blocked  <- function(mh_p_step = mh_p_step, mh_pprior_sd = mh_pprior_sd, 
                             p_1to0 = p1to0, p_0to1 = p_0to1, 
                             p_curr, occ_curr,
                             occur_prior_probs, probobs_curr, probobs_others, occur_others,
                             curr_inter, focus, 
                             detected){
  
  # Extract important structures
  p <- 2
  num_obs <- nrow(p_curr) # number of species of this type
  num_studies <- ncol(p_curr)
  pipj <- outer(probobs_curr, probobs_others)
  
  # Set up storage
  accepted = matrix(0, nrow = num_obs, ncol = num_studies)
  
    
  # Sample new occurrence and probability
  for(st in 1:num_studies){ # Loop over studies
      
      #### 1. Set up useful quantities for this study: all plants
      pi_prior_st <- occur_prior_probs[, st]
      p_curr_st  <- p_curr[,st]
      occur_curr_st <- occ_curr[,st]
      
      ### 2. Sample new p and new occurrence indicator for all nodes in this study, independently
      p_prop_st <- rtruncnorm(n = 1,  mean = p_curr_st, sd = mh_p_step, a = 0, b = 1) # Sample new p from trunc norm
      prop_probs_st <- ifelse(occur_curr_st == 1, 1-p_1to0, p_0to1)
      occur_prop_st <- rbinom(num_obs, 1, prop_probs_st)
      
      
      ### 3.  Accept or reject both p, occur for each node within study st
      
      # Set up useful values for both current and proposal values
      occ_others_st <- outer(rep(1,num_obs), occur_others[, st] )
      prod1 <- apply((1 - pipj) ^ (curr_inter * focus[, , st]* occ_others_st), 1, prod)   # Take product over other species
      
      ## Proposal values
      log_prob1_prop <- loglik(p_prop_st, prod1, occur_prop_st) # p(occ | pi, ...)
      log_prob_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))   # p(pi | ...) just the trunc normal prior
      log_prop_pi_prop <- log(dtruncnorm(x = p_prop_st, a = 0, b = 1, mean = p_curr_st, sd = mh_p_step)) # Asymmetric proposal for pi (trunc norm centered at p_prev)
      log_prop_occ_prop <- log(prop_probs_st^occur_prop_st * (1-prop_probs_st)^(1-occur_prop_st)) # Asymmetric proposal for occ (indicator switch)
      
      ## Current values
      log_prob1_curr <- loglik(p_curr_st, prod1, occur_curr_st) # p(occ | pi, ...)
      log_prob_pi_curr <- log(dtruncnorm(x = p_curr_st, a = 0, b = 1, mean = pi_prior_st, sd = mh_pprior_sd))   # p(pi | ...) just the trunc normal prior
      log_prop_pi_curr <- log(dtruncnorm(x = p_curr_st, a = 0, b = 1, mean = p_prop_st, sd = mh_p_step)) # Asymmetric proposal for pi (trunc norm centered at p_prev)
      prop_probs_st_curr <- ifelse(occur_prop_st == 1, 1-p_1to0, p_0to1) # Prob of moving from proposal to current
      log_prop_occ_curr <- log(prop_probs_st_curr^occur_curr_st * (1-prop_probs_st_curr)^(1-occur_curr_st)) # Asymmetric proposal for occ (indicator switch)
      
      ## Acceptance probability
      AP <- log_prob1_prop + log_prob_pi_prop + log_prop_pi_curr + log_prop_occ_curr # changed these two lines 
      AP <- AP - (log_prob1_curr + log_prob_pi_curr + log_prop_pi_prop + log_prop_occ_prop) 
      update <- log(runif(num_obs)) < AP
      
      ## Update stored quantities for the selected species
      p_curr[update, st] <- p_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      accepted[update, st] <- 1 
      
    } # End loop over studies st
    
    # Now, if a species was detected to interact in a study, set their occurrence
    # probability to 1
    # also set their occurrence to 1 (NEW - didn't do this previously)
    wh_detected <- which(detected == 1)
    p_curr[wh_detected] <- 1
    occ_curr[wh_detected] <- 1
    
  ret <- list(p_curr = p_curr, occ_curr = occ_curr, accepted = accepted)
  return(ret)
  
}

