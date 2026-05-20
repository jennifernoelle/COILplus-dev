#' Trait matching for species interactions
#'
#' Fits univariate regressions of the logit of posterior interaction
#' probability samples on covariates and uses permutation to obtain the
#' null distribution of no association. Returns both R-squared (unsigned
#' importance) and mean posterior correlation (signed direction) for each
#' covariate.
#'
#' @param B Integer. Number of permutation resamples. Default is 500.
#' @param mod_pL1s Array of posterior samples of fitted interaction
#'   probabilities, as returned by \code{MCMC}. Dimensions: [Nsims, nB, nP].
#' @param Xs List of posterior imputed covariate values for the first set of
#'   species, as returned by \code{MCMC}. Each element has rows equal to the
#'   number of posterior samples and columns equal to the number of missing
#'   observations for that covariate.
#' @param Ws Same as \code{Xs} but for the second set of species.
#' @param obs_X Design matrix of covariates for the first set of species
#'   (species x covariates).
#' @param obs_W Design matrix of covariates for the second set of species
#'   (species x covariates).
#' @param obs_only Logical. If \code{TRUE}, only species with fully observed
#'   covariates are used. Defaults to \code{FALSE}.
#' @param ncores Integer. Number of cores for parallel processing over
#'   posterior samples. Default is 1 (sequential). Values greater than 1
#'   use a \code{doParallel} cluster.
#'
#' @return A list with elements:
#' \describe{
#'   \item{rsq_resampling_X}{Permutation R-squared for covariates of species
#'     type 1. Matrix [B x covariates].}
#'   \item{rsq_resampling_W}{Permutation R-squared for covariates of species
#'     type 2. Matrix [B x covariates].}
#'   \item{rsq_obs_X}{Observed mean R-squared for each covariate of species
#'     type 1.}
#'   \item{rsq_obs_W}{Observed mean R-squared for each covariate of species
#'     type 2.}
#'   \item{corr_obs_X}{Posterior mean correlation between each species type 2
#'     and each covariate of species type 1. Matrix [nP x covariates].}
#'   \item{corr_obs_W}{Posterior mean correlation between each species type 1
#'     and each covariate of species type 2. Matrix [nB x covariates].}
#' }
#'
#' @export
#'
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doParallel registerDoParallel
TraitMatching <- function(B = 500, mod_pL1s, Xs, Ws, obs_X, obs_W,
                           obs_only = FALSE, ncores = 1) {

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

  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    `%op%` <- `%dopar%`
  } else {
    `%op%` <- `%do%`
  }

  cat('Covariates of first set of species.\n')

  for (mm in 1:sum_pB) {
    cat('\n Covariate', mm, '\n')

    res_list <- foreach(ss = 1:Nsims, .combine = rbind, .packages = "stats") %op% {
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

      rsq_obs <- mean(sapply(1:nP, function(jj) cor(this_response[, jj], this_cov)^2))

      rsq_perm <- numeric(B)
      for (bb in 1:B) {
        permuted_cov <- sample(this_cov)
        rsq_perm[bb] <- mean(sapply(1:nP, function(jj) cor(this_response[, jj], permuted_cov)^2))
      }

      c(rsq_obs, rsq_perm)
    }

    rsq_obs_X[, mm] <- res_list[, 1]
    rsq_resampling_X[, , mm] <- t(res_list[, -1])

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

    res_list <- foreach(ss = 1:Nsims, .combine = rbind, .packages = "stats") %op% {
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

      rsq_obs <- mean(sapply(1:nB, function(ii) cor(this_response[ii, ], this_cov)^2))

      rsq_perm <- numeric(B)
      for (bb in 1:B) {
        permuted_cov <- sample(this_cov)
        rsq_perm[bb] <- mean(sapply(1:nB, function(ii) cor(this_response[ii, ], permuted_cov)^2))
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
  names(dimnames(rsq_resampling_X)) <- c('boot', 'covariate')
  dimnames(rsq_resampling_W) <- list(boot = 1:B, cov = 1:sum_pP)
  names(dimnames(rsq_resampling_W)) <- c('boot', 'covariate')
  names(rsq_obs_X) <- 1:sum_pB
  names(rsq_obs_W) <- 1:sum_pP
  dimnames(corr_obs_X) <- list(species = 1:nP, cov = 1:sum_pB)
  names(dimnames(corr_obs_X)) <- c('species', 'covariate')
  dimnames(corr_obs_W) <- list(species = 1:nB, cov = 1:sum_pP)
  names(dimnames(corr_obs_W)) <- c('species', 'covariate')

  return(list(
    rsq_resampling_X = rsq_resampling_X,
    rsq_resampling_W = rsq_resampling_W,
    rsq_obs_X = rsq_obs_X,
    rsq_obs_W = rsq_obs_W,
    corr_obs_X = corr_obs_X,
    corr_obs_W = corr_obs_W
  ))
}
