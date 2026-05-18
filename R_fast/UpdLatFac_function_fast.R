#' Updating the model's latent factors
#' This version improves efficiency and accuracy by avoiding explicit inversions
#' 
#' MCMC sampling of latent factors for either set of species.
#' 
#' @param latfac Current values of the latent factors to be updated. Matrix
#' with number of rows equal to the number of species and number of columns
#' equal to the number of factors.
#' @param latfac_others Current values of the latent factors fot the other set
#' of species, not to be updated. Matrix. Rows correspond to species and
#' columns correspond to latent factors.
#' @param probobs Vector including the values for the probability of observing
#' a given species interaction within a study. Should be of length equal to the
#' number of species among this set.
#' @param coefs_probobs Coefficients of the latent factors for the normal model
#' on the logit of probability of observing a given species. If these
#' coefficients are all NA, then this implies that no bias correction for
#' geographical and taxonomical bias is performed, so this part of the model
#' does not inform the updates of the latent factors. Should be of length equal
#' to the number of latent factors used plus 1 for the intercept.
#' @param var_probobs Numeric. Variance in the model for the probability of
#' detecting an interaction. Used only if performing bias correction.
#' @param obs_covs  Matrix of observed traits. Rows correspond to species and
#' columns correspond to covariates. Continuous covariates should be included
#' first.
#' @param omega_obs_covs The Polya-Gamma omegas from the models of binary
#' covariates. Matrix. Rows correspond to units, columns to binary traits.
#' @param num_covs Vector of length 2 including the number of continuous and
#' binary traits.
#' @param coefs_covs Matrix of latent factor coefficients in the traits models.
#' Rows correspond to traits and columns to factors.
#' @param var_covs Vector of length equal to the number of continuous traits
#' including the residual variance of the corresponding model.
#' @param curr_inter Matrix with rows corresponding to the set of species whose
#' latent factors we are updating and columns corresponding to the other set of
#' species. Entries are 0 or 1.
#' @param coefs_inter Vector including the coefficients in the network model.
#' @param omega_inter The omegas generated from the Polya-Gamma for the model
#' of the true interactions L. The rows should correspond to the units whose
#' latent factors are updated using this function.
#' @param prior_S_inv The inverse of the current correlation matrix for the
#' latent factors that are being updated.
#' @param cut_feed Logical. Whether we should cut the feedback from the
#' interaction and detectability submodels into the latent factors. Defaults to
#' FALSE, in which case the full posterior is considered.
#' 
UpdLatFac_fast <- function(latfac, latfac_others,
                      probobs, coefs_probobs, var_probobs,
                      obs_covs, omega_obs_covs, num_covs, coefs_covs,
                      var_covs,
                      curr_inter, coefs_inter, omega_inter,
                      prior_S_inv, cut_feed = FALSE) {
  
  
  ret <- latfac
  Hval <- dim(latfac)[2]
  num_obs <- dim(latfac)[1]
  num_obs_other <- dim(latfac_others)[1]
  
  if (any(dim(omega_inter) != c(num_obs, num_obs_other))) {
    stop('Wrong dimensions of omega_inter.')
  }
  
  # If the coefficients for the observation model are NA, it implies that bias
  # correction is not performed.
  bias_cor <- !(all(is.na(coefs_probobs)))

  
  for (hh in 1 : Hval) {
    
    # Attaching the column of 1s to the latent factors:
    des_mat <- cbind(1, ret[, - hh])
    
    
    # Part 1(A): Combining the prior and the likelihood from continuous traits.
    # It leads to an updated normal distribution with parameters mpart, Spart.

  
    # Partial covariance matrix.
    coef_diag <- sum(coefs_covs[1 : num_covs[1], hh + 1] ^ 2 / var_covs)
    Spart_inv <- coef_diag * diag(num_obs) + prior_S_inv

    # Partial mean vector.
    mpart <- rep(0, num_obs)
    
    cont <- num_covs[1]
    bin <- num_covs[2]
    
    if (cont > 0) {
      # ensure matrix types (handles cont == 1)
      Xc     <- as.matrix(obs_covs[, seq_len(cont), drop = FALSE])                 # (n x cont)
      Cc_noh <- as.matrix(coefs_covs[seq_len(cont), -(hh + 1), drop = FALSE])      # (cont x H)
      P_cont <- des_mat %*% t(Cc_noh)                                              # (n x cont)
      resid_c <- Xc - P_cont                                                       # (n x cont)
      
      # weights: length cont numeric
      w_c <- as.numeric(coefs_covs[seq_len(cont), hh + 1] / var_covs[seq_len(cont)])
      
      # add weighted column-sum
      mpart <- mpart + drop(resid_c %*% w_c)
    }

    
    # Part 1(B): If bias correction is performed and the feedback is allowed, 
    # update the prior and continuous traits to include the model for the
    # probability of observing.
    
    if (bias_cor & !cut_feed) {
      
      coef_diag <- coefs_probobs[hh + 1] ^ 2 / var_probobs
      Spart_inv <- Spart_inv + coef_diag * diag(num_obs)
      
      pred <- des_mat %*% matrix(coefs_probobs[- (hh + 1)], ncol = 1)
      resid <- logit(probobs) - pred
      mpart <- mpart + coefs_probobs[hh + 1] * resid / var_probobs
      
    }
    
    
    # ----- Part 2: build posterior precision A = Spart_inv + diag(additions) -----
    # Updating the variance to include the Polya-Gamma components:
    # model for the true interactions L, and models for binary covariates.
    
    A <- Spart_inv
    diag_add <- numeric(num_obs)
    
    # Feedback (interaction) diagonal term: if feedback is allowed, we add the terms for interaction model
    if (!cut_feed) {
      v2 <- latfac_others[, hh]^2
      #row_add <- rowSums(sweep(omega_inter, 2, v2, `*`))
      row_add <- as.numeric(omega_inter %*% v2)  # same result, faster
      diag_add <- diag_add + (coefs_inter[hh + 1]^2) * row_add
    }
    
    # Binary covariate diagonal terms: adding the terms for the binary covariates
    if (bin > 0) {
      w_bin2  <- as.numeric(coefs_covs[cont + seq_len(bin), hh + 1]^2)  # length bin
      # adds: sum_j (use_coef_j^2 * omega[,j]) to the diagonal
      diag_add <- diag_add + as.vector(omega_obs_covs %*% w_bin2)
    }

    diag(A) <- diag(A) + diag_add
    
    # ----------------------- Part 3: construct b for posterior mean ----------------------#
    
    # Updating the mean of the conditional posterior to include the
    # Polya-Gamma components for the model for the true interactions L, and
    # the models for binary covariates.
    # This part is the same as old new_m assembly except we do not pre-multiply by covar
    
    # Starting from partial information
    b <- mpart
    
    # Adding the part corresponding to the L model
    if (!cut_feed) {
      coefs_excl <- coefs_inter[-c(1, hh + 1)]
      if (length(coefs_excl)) {
        A_noh <- ret[, -hh, drop = FALSE] *
          matrix(coefs_excl, nrow = num_obs, ncol = length(coefs_excl), byrow = TRUE)
        pred_noh <- A_noh %*% t(latfac_others[, -hh, drop = FALSE])
      } else {
        pred_noh <- matrix(0, nrow = num_obs, ncol = num_obs_other)
      }
      pred  <- coefs_inter[1] + pred_noh
      resid <- (curr_inter - 0.5) / omega_inter - pred
      step  <- (omega_inter * resid) %*% latfac_others[, hh, drop = FALSE]
      b <- b + as.numeric(coefs_inter[hh + 1] * step)
    }
    
    # Adding the part corresponding to the binary covariates #
    if (num_covs[2] > 0) {
      for (jj in 1:num_covs[2]) {
        use_coef <- coefs_covs[num_covs[1] + jj, hh + 1]
        pred <- des_mat %*% coefs_covs[num_covs[1] + jj, - (hh + 1)]
        resid <- (obs_covs[, num_covs[1] + jj] - 0.5) / omega_obs_covs[, jj] - pred
        b <- b + use_coef * resid * omega_obs_covs[, jj]
      }
    }
    
    
    # ----- Part 4: Solve instead of invert; sample without building covariance ----------#
    
    # Legacy version: explicitly form Sigma = A^{-1} and its Cholesky — slow and numerically unstable.
    # Current version: avoid forming Sigma by using a Cholesky of A and a few triangular solves.
    
    # Compute posterior mean mu via Cholesky and two triangular solves
    R  <- chol(A)                                    # A = R'R; chol() in R returns upper triangular R
    y  <- forwardsolve(t(R), b, upper.tri = FALSE)   # solve R' y = b
    mu <- backsolve(R, y, upper.tri = TRUE)          # solve R mu = y
    
    # Sample noise: draw z ~ N(0, I), then transform to have covariance A^{-1}
    z <- rnorm(num_obs)
    e <- backsolve(R, z, upper.tri = TRUE)           # e ~ N(0, A^{-1})
    
    ret[, hh] <- as.numeric(mu + e)
    
  }
  
  return(ret)
  
}


