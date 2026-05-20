#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector row_logprod_mask_idx_slice(const Rcpp::IntegerMatrix& curr_inter,
                                               const Rcpp::IntegerMatrix& focus_st,
                                               const Rcpp::IntegerVector& i_idx,
                                               const Rcpp::IntegerVector& j_idx,
                                               const Rcpp::NumericMatrix& log1m_pipj) {
  const int n = curr_inter.nrow();
  Rcpp::NumericVector out(n);  // default 0 for rows we won't touch
  
  const int ni = i_idx.size();
  const int nj = j_idx.size();
  
  for (int a = 0; a < ni; ++a) {
    int i = i_idx[a] - 1; // R->C++
    double acc = 0.0;
    for (int b = 0; b < nj; ++b) {
      int j = j_idx[b] - 1;
      if (curr_inter(i,j) && focus_st(i,j)) {
        acc += log1m_pipj(i,j);
      }
    }
    out[i] = acc; // log(product); untouched rows stay 0
  }
  return out;
}
