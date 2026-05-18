# This function includes block updates of occurrence indicators and probabilities
# Efficiency updates relative to legacy version

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
#' @param focus_rows_idx_list List of length `num_studies`. For study `st`,
#'   an integer vector of row indices `i` (1-based) such that the
#'   study-specific focus mask has at least one focused counterpart in that row:
#'   `any(focus[i, , st] != 0)`. This pre-filter limits rowwise work to only
#'   rows that could contribute. Indices must be unique, sorted, and within
#'   `1:nrow(curr_inter)`.
#'
#' @param focus_cols_idx_list List of length `num_studies`. For study `st`,
#'   anint eger vector of column indices `j` (1-based) with at least one
#'   focused counterpart: `any(focus[ , j, st] != 0)`. This pre-filter limits
#'   columnwise work to only columns that could contribute. Indices must be
#'   unique, sorted, and within `1:ncol(curr_inter)`.
#'
#' @param focus_int_list List of length `num_studies`. For study `st`, an
#'   integer (0/1) matrix with the same dimensions as `curr_inter`, equal to
#'   `focus[ , , st]` coerced to integer (i.e., `storage.mode(F) <- "integer"`).
#'   This is the per-study focus mask used by the C++ routine
#'   `row_logprod_mask_idx_slice()`. Dimensions must exactly match
#'   `nrow(curr_inter)` × `ncol(curr_inter)`.
#'
#' @section Details on orientation:
#'   The row/column meaning is with respect to `curr_inter` passed to the
#'   function. For the P-block (plants-as-rows), you typically use
#'   `curr_inter = t(this_L)` and compute the lists from
#'   `focus_perm <- aperm(focus, c(2, 1, 3))` so that rows/cols line up.
#'   For the B-block (vertebrates-as-rows), use `curr_inter = this_L` and the
#'   original `focus`.
#'
#'   Invariants expected by the function:
#'   - `length(focus_rows_idx_list) == length(focus_cols_idx_list) ==
#'     length(focus_int_list) == num_studies`
#'   - For each `st`, `nrow(focus_int_list[[st]]) == nrow(curr_inter)` and
#'     `ncol(focus_int_list[[st]]) == ncol(curr_inter)`.
#'   - No `NA` in the index vectors; indices are 1-based.
#'
#' @note Precomputing these lists once per run (or when `focus` changes)
#'   avoids repeated scans of the 3D `focus` array inside the MCMC loop and
#'   enables the C++ masked rowwise product to work on minimal slices.


UpdOccurP_blocked_fast <- function(
    mh_p_step, mh_pprior_sd, p_1to0, p_0to1,
    p_curr, occ_curr,
    occur_prior_probs, probobs_curr, probobs_others, occur_others,
    curr_inter, 
    focus_rows_idx_list, focus_cols_idx_list, focus_int_list,
    logZ_prior_mat = NULL, # precomputed prior
    detected   # matrix [num_obs x num_studies] (pass detected_P / detected_B)
    ){
  
  
  # sizes
  num_obs     <- nrow(p_curr)
  num_studies <- ncol(p_curr)
  
  # normalize detected to logical matrix (per-study), preserving original semantics
  if (!(is.matrix(detected) && all(dim(detected) == c(num_obs, num_studies))))
    stop("`detected` must be a matrix with dim = nrow(p_curr) x ncol(p_curr)")
  detected_mat <- detected != 0L
  
  # per-call things (depend only on pi/pj)
  log1m_pipj <- log1p(-outer(probobs_curr, probobs_others))  # n_i x n_j
  
  accepted <- matrix(0L, nrow = num_obs, ncol = num_studies)
  
  
  for (st in seq_len(num_studies)) {

    
    ## ACTIVE MASK (per-study): rows we may update
    active_idx <- which(!detected_mat[, st])
    
    ## --- propose (new p, new occ) ---
    pi_prior_st    <- occur_prior_probs[, st]
    p_curr_st      <- p_curr[, st]
    occur_curr_st  <- occ_curr[, st]
    
    # draw for ALL rows (same as old behavior)
    p_prop_full <- rtruncnorm(n = 1, mean = p_curr_st, sd = mh_p_step, a = 0, b = 1)
    occ_prop_full <- rbinom(
      n = num_obs, size = 1L,
      prob = ifelse(occur_curr_st == 1L, 1 - p_1to0, p_0to1)
    )
    
    # then mask: inactive rows keep current values
    p_prop_st     <- p_curr_st
    occur_prop_st <- occur_curr_st
    if (length(active_idx)) {
      p_prop_st[active_idx]     <- p_prop_full[active_idx]
      occur_prop_st[active_idx] <- occ_prop_full[active_idx]
    }

    
    ## --- prod1 (restrict rows to active) ---

    occ_idx <- intersect(which(occur_others[, st] == 1L), focus_cols_idx_list[[st]])
    if (length(occ_idx) == 0L) {
      log_prod1 <- numeric(num_obs)  # zeros
    } else {
      i_idx    <- intersect(focus_rows_idx_list[[st]], active_idx)
      if (length(i_idx)) {
        log_prod1 <- row_logprod_mask_idx_slice(
          curr_inter, focus_int_list[[st]],
          as.integer(i_idx), as.integer(occ_idx), log1m_pipj
        )
      } else {
        log_prod1 <- numeric(num_obs)
      }
    }
    
    ## --- log-likelihood & priors (only active rows) ---

    # logZ_prior     <- log(pnorm((1 - pi_prior_st)/mh_pprior_sd) - pnorm((-pi_prior_st)/mh_pprior_sd))
    logZ_prior <- if (is.null(logZ_prior_mat)) {
      log(pnorm((1 - pi_prior_st)/mh_pprior_sd) - pnorm((-pi_prior_st)/mh_pprior_sd))
    } else {
      logZ_prior_mat[, st]
    }
    logZ_step_curr <- log(pnorm((1 - p_curr_st)/mh_p_step)      - pnorm((-p_curr_st)/mh_p_step))
    logZ_step_prop <- log(pnorm((1 - p_prop_st)/mh_p_step)      - pnorm((-p_prop_st)/mh_p_step))
    
    APa <- NULL  # acceptance stat for active rows only
    if (length(active_idx)) {
      # (keep your existing lp_prop/lp_curr + truncated-normal + Bernoulli pieces)
      lp_prop <- fast_loglik_log(p_prop_st[active_idx], log_prod1[active_idx], occur_prop_st[active_idx])
      lp_curr <- fast_loglik_log(p_curr_st[active_idx],  log_prod1[active_idx], occur_curr_st[active_idx])
      
      log_prob_pi_prop <- log_dtrunc01(p_prop_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd, logZ_prior[active_idx])
      log_prob_pi_curr <- log_dtrunc01(p_curr_st[active_idx],  pi_prior_st[active_idx], mh_pprior_sd, logZ_prior[active_idx])
      log_prop_pi_prop <- log_dtrunc01(p_prop_st[active_idx],  p_curr_st[active_idx],   mh_p_step,    logZ_step_curr[active_idx])
      log_prop_pi_curr <- log_dtrunc01(p_curr_st[active_idx],  p_prop_st[active_idx],   mh_p_step,    logZ_step_prop[active_idx])
      
      eps <- 1e-12
      pp  <- p_0to1 + occur_curr_st[active_idx] * ((1 - p_1to0) - p_0to1)
      ppc <- p_0to1 + occur_prop_st[active_idx] * ((1 - p_1to0) - p_0to1)
      pp  <- pmin(pmax(pp,  eps), 1 - eps)
      ppc <- pmin(pmax(ppc, eps), 1 - eps)
      log_prop_occ_prop <- occur_prop_st[active_idx]*log(pp)  + (1-occur_prop_st[active_idx])*log1p(-pp)
      log_prop_occ_curr <- occur_curr_st[active_idx]*log(ppc) + (1-occur_curr_st[active_idx])*log1p(-ppc)
      
      # Diagnostics printing
      if (anyNA(p_prop_st[active_idx]) || anyNA(p_curr_st[active_idx])) {
        warning(sprintf("[UpdOccurP] NA in p_* (st=%d): prop=%d curr=%d",
                        st, sum(is.na(p_prop_st[active_idx])), sum(is.na(p_curr_st[active_idx]))))
      }
      if (any(!is.finite(log_prod1[active_idx]))) {
        warning(sprintf("[UpdOccurP] non-finite log_prod1 (st=%d): %d",
                        st, sum(!is.finite(log_prod1[active_idx]))))
      }
      if (!(mh_p_step > 0 && mh_pprior_sd > 0)) stop("mh_p_step and mh_pprior_sd must be > 0")
      
      
      APa <- (lp_prop + log_prob_pi_prop + log_prop_pi_curr + log_prop_occ_curr) -
        (lp_curr + log_prob_pi_curr + log_prop_pi_prop + log_prop_occ_prop)
      # Treat any non-finite AP as automatic reject
      APa[!is.finite(APa)] <- -Inf
      
    }
    
    ## --- accept/reject + writeback ---
    # NOTE: `accepted` marks MH acceptances only on non-detected rows (detected==1 are skipped).
    # This is intentional for speed; detected rows would have been enforced to 1's anyway.
    
    u_full <- runif(num_obs)                 # draw for ALL rows (preserves RNG stream)
    update <- rep(FALSE, num_obs)
    if (length(active_idx)) {
      update[active_idx] <- is.finite(APa) & (log(u_full[active_idx]) < APa)
    }
    
    if (any(update)) {
      p_curr[update, st]   <- p_prop_st[update]
      occ_curr[update, st] <- occur_prop_st[update]
      accepted[update, st] <- 1L
    }
  }
  
  # no row-wise “enforce all studies” block here (preserves original sparsity)
  
  list(
    p_curr   = p_curr,
    occ_curr = occ_curr,
    accepted = accepted
  )
}
