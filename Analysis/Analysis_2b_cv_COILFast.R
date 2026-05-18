# Performs 10-fold CV with COIL

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo

# Save results using convention: res_date_i.rda
date <- 'COIL_exp'

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path_base <- 'Results/'
# Where the functions are available:
source_path <- 'R_fast/'

# Create the results folder
ifelse(!dir.exists(file.path(save_path_base, date)), dir.create(file.path(save_path_base, date)), FALSE)
save_path <- paste0('Results/', date, '/')


# ---------------------- STEP 0: Load functions, packages, data -------------- #

# Unchanged legacy functions
source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))

# Updated fast functions
Rcpp::sourceCpp(paste0(source_path, "row_logprod_mask_idx_slice.cpp"))
source(paste0(source_path, 'UpdLatFac_function_fast.R'))
source(paste0(source_path, 'UpdOccurP_function_blocked_fast.R'))
source(paste0(source_path, 'Utils_OccurP.R'))
source(paste0(source_path, 'UpdOccur_function_fast.R'))
source(paste0(source_path, 'MCMC_fast.R'))

library(parallel)
library(foreach)
library(abind)
library(magrittr)
library(truncnorm)
# Loading the data:
load(paste0(data_path, 'Cu_phylo_sorted.dat'))
load(paste0(data_path, 'Cv_phylo_sorted.dat'))
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'F_obs.dat'))
load(paste0(data_path, 'Obs_X_sorted.dat')) # vertebrate traits
load(paste0(data_path, 'Obs_W_sorted.dat')) # plant traits
load(paste0(data_path, 'OP_full.dat')) # site level obs plants
load(paste0(data_path, 'OV_full.dat')) # site level obs vertebrates

## Rename for convenience
Cu <- Cu_sorted
Cv <- Cv_sorted

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]

# -------------- STEP 1: Specifications. ------------ #

repetitions <- 10 # Number of cv replications
n.cv <- 100 # Number of heldout observations

Nsims <- 5 # Recommended 500 
thin <-  1 # Recommended 20 
burn <-  5 # Recommended 10000 
use_H <- 10 
theta_inf <- 0.01

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

# sampling <- NULL
sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                 gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                 rV = TRUE, miss_X = TRUE, miss_W = TRUE, O_V = TRUE,
                 O_P = TRUE, p_OV = FALSE, p_OP = FALSE)

start_values <- NULL
block_sampleOccP <- FALSE
bias_cor <- TRUE 


# --------------------------- STEP 2: CROSS VALIDATION ------------------------- #


# We recommend running in parallel 
  
for(rr in 1:repetitions){
  set.seed(rr)
  cat("\n Running cross validation fold ", rr)
  
  # Matrix that chooses n.cv recorded interactions to be held-out.
  set_out <- matrix(0, nrow = nV, ncol = nP)
  set_out[sample(which(comb_A == 1), n.cv)] <- 1
  
  
  # Getting the indices that were zero-ed out.
  cv_indices <- matrix(NA, nrow = n.cv, ncol = 2)
  wh <- which(set_out == 1)
  for (ii in 1 : n.cv) {
    row_ii <- wh[ii] %% nV
    row_ii <- ifelse(row_ii == 0, nV, row_ii)
    col_ii <- ceiling(wh[ii] / nV)
    cv_indices[ii, ] <- c(row_ii, col_ii)
  }
  
  # Zero-ing out corresponding entries in A. 
  use_A <- obs_A
  for (ii in 1 : n.cv) {
    use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Also set Focus to zero for these observations?
  use_F <- obs_F
  for (ii in 1 : n.cv) {
    use_F[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  mcmc <- MCMC(obs_A = use_A, focus = use_F, p_occur_V = O_V, p_occur_P = O_P, obs_X, obs_W, Cu, Cv,
               Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
               bias_cor = TRUE, theta_inf = 0.01,
               mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
               mh_pprior_sd = 0.1, mh_p_step = 0.1,
               stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
               prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
               prior_sigmasq = c(1, 1), start_values = NULL,
               sampling = sampling,
               cut_feed = FALSE,
               block_sampleOccP = block_sampleOccP)
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean)
  save(pred, file = paste0(save_path, 'pred_', date, '_', rr, '.dat'))
  rm(pred)
  
  save(cv_indices, file = paste0(save_path, 'cv_indices_', date, '_', rr, '.dat'))
  rm(cv_indices)
  
}