#' Function that updates the indicator of occurence. Corrects errors in original version
#' 
#' Each species might occur in the study area or not. Even though probabilities
#' of occurence can be provided, the detected interactions also inform us of
#' whether the species occur.
#' 
#' @param detected Matrix. Rows are species, columns are studies. Indicator of
#' whether the species were detected to be present in the given study. Derived
#' from the provided occurrence matrix, and checked against obs_A. Note that
#' while most studies only report interacting species, some studies do report
#' present but not interacting species. This often occurs when subsetting. 
#' @param occur Matrix. Rows are species, columns are studies. It includes the
#' prior probabilities of occurence.
#' @param occur_others. Matrix Rows are species, columns are studies. It includes
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
#' @param focus_rows_idx_list List of length `num_studies`. For study `st`,
#'   an integer vector of row indices `i` (1-based) such that the
#'   study-specific focus mask has at least one focused counterpart in that row:
#'   `any(focus[i, , st] != 0)`. This pre-filter limits rowwise work to only
#'   rows that could contribute. Indices must be unique, sorted, and within
#'   `1:nrow(curr_inter)`.
#'
#' @param focus_cols_idx_list List of length `num_studies`. For study `st`,
#'   an integer vector of column indices `j` (1-based) with at least one
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


# Fast occurrence-only updater (no p updates)
# deps: row_logprod_mask_idx_slice (your C++ helper)
UpdOccur_fast <- function(detected, occur, occur_others,
                          probobs_curr, probobs_others,
                          curr_inter,
                          focus_rows_idx_list, focus_cols_idx_list, focus_int_list,
                          log1m_pipj = NULL) {
  num_obs     <- nrow(occur)
  num_studies <- ncol(occur)
  if (ncol(occur_others) != num_studies) stop("occur_others ncol mismatch")
  if (!is.matrix(detected) || any(dim(detected) != dim(occur)))
    stop("detected must be a matrix with same dim as occur")
  
  # per-call: log(1 - p_i p_j)
  if (is.null(log1m_pipj)) {
    log1m_pipj <- log1p(-outer(probobs_curr, probobs_others))  # n_i x n_j
  }
  
  new_O <- matrix(0L, num_obs, num_studies)
  eps <- 1e-12
  
  for (st in seq_len(num_studies)) {
    # columns that matter in this study (occur==1 AND in focus)
    occ_cols <- focus_cols_idx_list[[st]]
    occ_idx  <- if (length(occ_cols)) occ_cols[ occur_others[occ_cols, st] == 1L ] else integer(0)
    
    # rows to compute (only rows with any focus this study)
    if (length(occ_idx)) {
      i_idx <- focus_rows_idx_list[[st]]
      if (length(i_idx)) {
        log_prod1 <- row_logprod_mask_idx_slice(
          curr_inter,
          focus_int_list[[st]],       # integer matrix (0/1)
          as.integer(i_idx),
          as.integer(occ_idx),
          log1m_pipj
        )
      } else {
        log_prod1 <- numeric(num_obs)
      }
    } else {
      log_prod1 <- numeric(num_obs)
    }
    
    # prior odds + likelihood in logit space: eta = logit(occur) + log_prod1
    occ_prior <- pmin(pmax(occur[, st], eps), 1 - eps)
    eta <- log(occ_prior / (1 - occ_prior)) + log_prod1
    
    # stable sigmoid
    s <- log1p(exp(-abs(eta))) + pmax(eta, 0)
    prob1 <- exp(eta - s)
    
    # draw for ALL rows (preserve RNG stream vs old function)
    new_O[, st] <- rbinom(num_obs, 1L, prob1)
  }
  
  # enforce detected = 1
  new_O[detected == 1L] <- 1L
  new_O
}
