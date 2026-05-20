#' Trait matching for species interactions (parallelized): adds signed correlation
#' Fitting univariate regressions of the logit of fitted interaction
#' probability values on covariates and using permutation to get their
#' distribution under the null of no association.
#' 
#' @param B Numeric. Number of times to perform the resampling. Default is 500.
#' @param mod_pL1s Posterior samples for the fitted probabilities of the
#' interaction model. Returned by the function MCMC.
#' @param Xs Posterior samples of imputed values for the covariates of the
#' first set of species. Returned by the function MCMC. Should be a list of 
#' length equal to the number of covariates, each element should have rows equal
#' to the number of posterior samples across all chains, and columns equal to 
#' the number of missing observations for that covariate.
#' @param Ws Posterior samples of imputed values for the covariates of the
#' second set of species. Returned by the function MCMC.
#' @param obs_X The design matrix of covariates for the first set of species.
#' Number of rows is number of species and number of columns is number of
#' covariates.
#' @param obs_W The design matrix of covariates for the second set of species.
#' Number of rows is number of species and number of columns is number of
#' covariates.
#' @param obs_only Logical. If set to TRUE only the observations with observed
#' covariate values will be used. Defaults to FALSE.
#' 
#' @export
#' 
#' 

library(doParallel)
library(foreach)

TraitMatching_parallel <- function(B = 500, mod_pL1s, Xs, Ws, obs_X, obs_W, obs_only = FALSE, 
                                    ncores) {
  #ncores <- parallel::detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  on.exit(stopCluster(cl))  # Ensure cluster shuts down on exit/error
  
  Nsims <- dim(mod_pL1s)[1]
  nB <- dim(mod_pL1s)[2]
  nP <- dim(mod_pL1s)[3]
  sum_pB <- ncol(obs_X)
  sum_pP <- ncol(obs_W)
  
  logit_mod_pL1s <- logit(mod_pL1s)
  
  rsq_resampling_X <- array(NA, dim = c(B, Nsims, sum_pB))
  rsq_resampling_W <- array(NA, dim = c(B, Nsims, sum_pP))
  rsq_obs_X <- array(NA, dim = c(Nsims, sum_pB))
  corr_obs_X <- array(NA, dim = c(nP, sum_pB))
  rsq_obs_W <- array(NA, dim = c(Nsims, sum_pP))
  corr_obs_W <- array(NA, dim = c(nB, sum_pP))
  
  cat('Covariates of first set of species.\n')
  
  for (mm in 1:sum_pB) {
    cat('\n Covariate', mm, '\n')
    
    res_list <- foreach(ss = 1:Nsims, .combine = rbind, .packages = "stats") %dopar% {
      this_response <- logit_mod_pL1s[ss, , ]
      this_cov <- obs_X[, mm]
      wh_na <- which(is.na(this_cov))
      
      if (length(wh_na) > 0) {
        if (obs_only) {
          this_cov <- this_cov[-wh_na]
          this_response <- this_response[-wh_na, ]
        } else {
          this_cov[wh_na] <- Xs[[mm]][ss, ]
        }
      }
      
      rsq <- sapply(1:nP, function(jj) cor(this_response[, jj], this_cov)^2)
      rsq_obs <- mean(rsq)
      
      rsq_perm <- numeric(B)
      for (bb in 1:B) {
        permuted_cov <- sample(this_cov)
        rsq_bb <- sapply(1:nP, function(jj) cor(this_response[, jj], permuted_cov)^2)
        rsq_perm[bb] <- mean(rsq_bb)
      }
      
      c(rsq_obs, rsq_perm)
    }
    
    rsq_obs_X[, mm] <- res_list[, 1]
    rsq_resampling_X[, , mm] <- t(res_list[, -1])
    
    # Correlations (not parallelized—fast)
    for (jj in 1:nP) {
      this_cov <- obs_X[, mm]
      wh_na <- which(is.na(this_cov))
      this_response_jj <- logit_mod_pL1s[, , jj]
      corr <- rep(NA, Nsims)
      
      for (ss in 1:Nsims) {
        if (length(wh_na) > 0) {
          if (obs_only) {
            this_cov <- this_cov[-wh_na]
            this_response_jj <- this_response_jj[, -wh_na]
          } else {
            this_cov[wh_na] <- Xs[[mm]][ss, ]
          }
        }
        corr[ss] <- cor(this_response_jj[ss, ], this_cov)
      }
      corr_obs_X[jj, mm] <- mean(corr)
    }
  }
  
  cat('Covariates of second set of species.\n')
  
  for (ll in 1:sum_pP) {
    cat('Covariate', ll, '\n')
    
    res_list <- foreach(ss = 1:Nsims, .combine = rbind, .packages = "stats") %dopar% {
      this_cov <- obs_W[, ll]
      wh_na <- which(is.na(this_cov))
      this_response <- logit_mod_pL1s[ss, , ]
      
      if (length(wh_na) > 0) {
        if (obs_only) {
          this_cov <- this_cov[-wh_na]
          this_response <- this_response[, -wh_na]
        } else {
          this_cov[wh_na] <- Ws[[ll]][ss, ]
        }
      }
      
      rsq <- sapply(1:nB, function(ii) cor(this_response[ii, ], this_cov)^2)
      rsq_obs <- mean(rsq)
      
      rsq_perm <- numeric(B)
      for (bb in 1:B) {
        permuted_cov <- sample(this_cov)
        rsq_bb <- sapply(1:nB, function(ii) cor(this_response[ii, ], permuted_cov)^2)
        rsq_perm[bb] <- mean(rsq_bb)
      }
      
      c(rsq_obs, rsq_perm)
    }
    
    rsq_obs_W[, ll] <- res_list[, 1]
    rsq_resampling_W[, , ll] <- t(res_list[, -1])
    
    for (jj in 1:nB) {
      this_cov <- obs_W[, ll]
      wh_na <- which(is.na(this_cov))
      this_response_jj <- logit_mod_pL1s[, jj, ]
      corr <- rep(NA, Nsims)
      
      for (ss in 1:Nsims) {
        if (length(wh_na) > 0) {
          if (obs_only) {
            this_cov <- this_cov[-wh_na]
            this_response_jj <- this_response_jj[, -wh_na]
          } else {
            this_cov[wh_na] <- Ws[[ll]][ss, ]
          }
        }
        corr[ss] <- cor(this_response_jj[ss, ], this_cov)
      }
      corr_obs_W[jj, ll] <- mean(corr)
    }
  }
  
  rsq_resampling_X <- apply(rsq_resampling_X, c(1, 3), mean)
  rsq_resampling_W <- apply(rsq_resampling_W, c(1, 3), mean)
  rsq_obs_X <- apply(rsq_obs_X, 2, mean)
  rsq_obs_W <- apply(rsq_obs_W, 2, mean)
  
  dimnames(rsq_resampling_X) <- list(boot = 1:B, cov = 1:sum_pB)
  dimnames(rsq_resampling_W) <- list(boot = 1:B, cov = 1:sum_pP)
  names(rsq_obs_X) <- 1:sum_pB
  names(rsq_obs_W) <- 1:sum_pP
  dimnames(corr_obs_X) <- list(species = 1:nP, cov = 1:sum_pB)
  dimnames(corr_obs_W) <- list(species = 1:nB, cov = 1:sum_pP)
  
  return(list(
    rsq_resampling_X = rsq_resampling_X,
    rsq_resampling_W = rsq_resampling_W,
    rsq_obs_X = rsq_obs_X,
    rsq_obs_W = rsq_obs_W,
    corr_obs_X = corr_obs_X,
    corr_obs_W = corr_obs_W
  ))
}
