# This file runs COIL across four chains using: 
# Expert-defined occurrence probabilities


# ---------------------------------- TO DO ------------------------------------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
# But you may need to supply it directly if you're running this in via slurm etc
# wd_path <- "YOUR PATH"
# setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'COIL_exp'

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path_base <- 'Results/' # Crete this directory
# Where the functions are available:
source_path <- 'R_fast/'

# Create the results folder
ifelse(!dir.exists(file.path(save_path_base, date)), dir.create(file.path(save_path_base, date)), FALSE)
save_path <- paste0(save_path_base, date, '/')

# ------ STEP 0: Some functions. --------- #


# Unchanged functions
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
source(paste0(source_path, 'Utils_OccurP.R')) # added fast_loglik
source(paste0(source_path, 'UpdOccur_function_fast.R'))
source(paste0(source_path, 'MCMC_fast.R'))



library(parallel)
library(foreach)
library(abind)
library(magrittr)
library(truncnorm)
library(BayesLogit)
library(mvnfast)

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

# Getting the combined network: interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]


# ---------------------------- STEP 1: Specifications. ----------------------- #

# Recommended and trial MCMC parameters
Nsims <- 10 # Recommended 500
thin <-  1 # Recommended 20 
burn <-  5 #Recommended 10000
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

sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                 gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                 rV = TRUE, miss_X = TRUE, miss_W = TRUE, 
                 O_V = TRUE, O_P = TRUE, # sample occurrence indicators
                 p_OV = FALSE, p_OP = FALSE # sample occurrence probabilities
                 )

start_values <- NULL
bias_cor <- TRUE # Performing bias correction


#------------------------------ FIT THE MODEL ---------------------------------#

# We recommend running four chains in parallel
n.chains <- 4 # Recommended: 4
for(cc in 1:n.chains){
  set.seed(cc)
  
  cat("\n Fitting chain ", cc)
  
  ## Create a unique filename for each interation of the parallel loop
  each_filename <- paste0('res_', date, '_', as.character(cc), '.dat')
  each_filepath <- file.path(save_path, each_filename)
  
  mcmc <- MCMC(obs_A, focus = obs_F, p_occur_V = O_V, p_occur_P = O_P, obs_X, obs_W, Cu, Cv,
               Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
               bias_cor = TRUE, theta_inf = 0.01,
               mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
               mh_pprior_sd = 0.1, mh_p_step = 0.1,
               stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
               prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
               prior_sigmasq = c(1, 1), start_values = NULL,
               sampling = sampling, cut_feed = FALSE, 
               save_logL = TRUE
               ) 
  
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  # Running mean of detection probabilities
  p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
  
  # Running mean of latent factors
  factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
  
  # Imputed values of missing covariates
  Xs <- mcmc$Xs
  Ws <- mcmc$Ws
  
  # Occurrence indicators, probabilities, acceptance rates (last two not relevant if we don't sample P)
  occ_plants <- list(OP_mean = mcmc$OP_mean, p_OPs = mcmc$p_OP_mean, p_accept = mcmc$p_OP_accepted) 
  occ_verts <- list(OV_mean = mcmc$OV_mean, p_VPs = mcmc$p_OV_mean, p_accept = mcmc$p_OV_accepted) 
  
  # Log likelihood
  logL <- mcmc$logL
  
  # Combining the results we are interested in to a list and saving:
  res <- list(all_pred = all_pred, logL = logL, correlations = correlations, 
              p_detect = p_detect, factors = factors, occ_plants = occ_plants, 
              occ_verts = occ_verts,  Xs = Xs, Ws = Ws)
  
  save(res, file = each_filepath)
  
  rm(res)
}


