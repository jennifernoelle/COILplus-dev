# Estimating variable importance based on interaction model

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# Set the number of cores avaialble and if you want to use the parallelized version
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
use_parallel <- TRUE


data_path <- 'ProcessedData/'
save_plots_path <- 'Results/'
source_path <- 'R_fast/'


# Save results using convention: res_date_i.rda
date <- "COILp_exp"
results_path <- paste0('Results/', date , '/')
save_path <- paste0('Results/', date, '/')

# ------ STEP 0: Some functions. --------- #

library(tidyverse)
library(grid)
library(gridExtra)

source(paste0(source_path, 'TraitMatching_function.R'))
source(paste0(source_path, 'TraitMatching_parallel_function.R'))
source(paste0(source_path, 'useful_functions.R'))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Loading the data:
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'Obs_X_sorted.dat')) # vertebrate traits
load(paste0(data_path, 'Obs_W_sorted.dat')) # plant traits


# Check for missing values in the covariates:
any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)

# Getting the sample sizes:
nV <- nrow(obs_X)
nP <- nrow(obs_W)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1


# --------------- STEP 1: Efficiently loading posterior samples ---------------#

# Number of MCMC chains
nchains <- 4

# Load all chains
all_pred <- vector("list", nchains)
all_X <- vector("list", nchains)
all_W <- vector("list", nchains)

for (ii in 1:nchains) {
  cat("\nChain =", ii)
  res <- loadRData(paste0(results_path, 'res_', date, '_', ii, '.dat')) 
  all_pred[[ii]] <- res$all_pred
  if (any_X_miss) all_X[[ii]] <- res$Xs
  if (any_W_miss) all_W[[ii]] <- res$Ws
}

# Set number of posterior samples and trait dimensions
Nsims <- dim(all_pred[[1]])[1]
pV <- if (any_X_miss) length(all_X[[1]]) else ncol(obs_X)
pP <- if (any_W_miss) length(all_W[[1]]) else ncol(obs_W)

# Combine predictions across chains
mod_pL1s <- array(NA, dim = c(nchains * Nsims, nV, nP))
for (cc in 1:nchains) {
  wh_entries <- Nsims * (cc - 1) + 1:Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}

# Construct imputed trait matrices if needed
Xs <- NULL
if (any_X_miss) {
  Xs <- vector("list", pV)
  for (p in 1:pV) {
    nmiss.p <- length(all_X[[1]][[p]])
    Xs[[p]] <- matrix(NA, nchains * Nsims, nmiss.p)
    for (cc in 1:nchains) {
      wh_entries <- Nsims * (cc - 1) + 1:Nsims
      Xs[[p]][wh_entries, ] <- all_X[[cc]][[p]]
    }
  }
}

Ws <- NULL
if (any_W_miss) {
  Ws <- vector("list", pP)
  for (p in 1:pP) {
    nmiss.p <- length(all_W[[1]][[p]])
    Ws[[p]] <- matrix(NA, nchains * Nsims, nmiss.p)
    for (cc in 1:nchains) {
      wh_entries <- Nsims * (cc - 1) + 1:Nsims
      Ws[[p]][wh_entries, ] <- all_W[[cc]][[p]]
    }
  }
}

# --------------- STEP 2: Performing trait matching ----------------- #

B <- 10 # 500 permutations recommended

# Dealing with extreme values for which logit(x) is infinite.
mod_pL1s[mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}

# Optional: reduce memory use by keeping fewer posterior draws
keep_idx <- seq(1, dim(mod_pL1s)[1], by = 1)  # Change by argument for thinning
mod_pL1s <- mod_pL1s[keep_idx, , ]


t1 <- Sys.time()
if(use_parallel){
  trait_match <- TraitMatching_parallel(B = B, mod_pL1s = mod_pL1s,
                                         Xs = NULL, Ws = NULL,  # Imputed values not used.
                                         obs_X = obs_X, obs_W = obs_W, obs_only = TRUE, 
                                         ncores = ncores)
}else{
  trait_match <- TraitMatching(B = B, mod_pL1s = mod_pL1s,
                                Xs = NULL, Ws = NULL,  # Imputed values not used.
                                obs_X = obs_X, obs_W = obs_W, obs_only = TRUE)
  
}


print(Sys.time() - t1)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W
corr_obs_X <- trait_match$corr_obs_X
corr_obs_W <- trait_match$corr_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X_', date, '.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W_', date, '.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W_', date, '.dat'))
save(corr_obs_X, file = paste0(save_path, 'corr_obs_X_', date, '.dat'))
save(corr_obs_W, file = paste0(save_path, 'corr_obs_W_', date, '.dat'))
