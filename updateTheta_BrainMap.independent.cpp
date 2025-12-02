// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List UpdateTheta_BrainMap_independent_cpp(
    const arma::mat& prior_mean,     // nV x nL
    const arma::mat& prior_var,      // nV x nL
    const arma::mat& BOLD,           // nV x nT
    const List& theta,               // list with A (nT x nL), nu0_sq (scalar)
    const arma::vec& C_diag,         // length nT
    Nullable<arma::mat> H = R_NilValue,
    Nullable<arma::mat> Hinv = R_NilValue,
    bool update_nu0sq = true,
    bool return_MAP = false,
    bool verbose = true
) {
  // dims
  int nV = BOLD.n_rows;
  int nT = BOLD.n_cols;
  int nL = prior_mean.n_cols;
  
  // extract theta$A and theta$nu0_sq
  arma::mat A = as<arma::mat>(theta["A"]); // nT x nL
  double nu0_sq = as<double>(theta["nu0_sq"]);
  
  // initialize outputs and accumulators
  arma::mat A_part1(nT, nL, fill::zeros);
  arma::mat A_part2(nL, nL, fill::zeros);
  
  arma::mat miu_s(nV, nL, fill::zeros);
  arma::cube miu_ssT(nL, nL, nV, fill::zeros); // slice = miu_ssT_v
  arma::mat var_s(nV, nL, fill::zeros);
  
  // precompute nu0C_inv, At_nu0Cinv, At_nu0Cinv_A
  arma::mat nu0C_inv = diagmat(1.0 / (C_diag * nu0_sq)); // nT x nT
  arma::mat At_nu0Cinv = A.t() * nu0C_inv;                // nL x nT
  arma::mat At_nu0Cinv_A = At_nu0Cinv * A;                // nL x nL
  
  // placeholders for last E_v_inv and Sigma_s_v (to match R's Estep which holds final values)
  arma::mat last_E_v_inv(nL, nL, fill::zeros);
  arma::mat last_Sigma_s_v(nL, nL, fill::zeros);
  
  // main loop over voxels
  for (int vv = 0; vv < nV; ++vv) {
    // y_v: row, s0_v: row, var row
    arma::rowvec y_v = BOLD.row(vv);           // 1 x nT
    arma::rowvec s0_v = prior_mean.row(vv);    // 1 x nL
    arma::rowvec var_v = prior_var.row(vv);    // 1 x nL
    
    // E-step
    arma::mat E_v_inv = diagmat(1.0 / var_v.t());           // nL x nL
    arma::mat Sigma_s_v = inv(E_v_inv + At_nu0Cinv_A);      // nL x nL
    arma::colvec miu_s_v = Sigma_s_v * 
      (At_nu0Cinv * y_v.t() + E_v_inv * s0_v.t());         // nL x 1
    arma::mat miu_ssT_v = miu_s_v * miu_s_v.t() + Sigma_s_v;// nL x nL
    
    // save posterior moments (for M-step of nu0_sq)
    miu_s.row(vv) = miu_s_v.t();
    for (int i = 0; i < nL; ++i) var_s(vv, i) = Sigma_s_v(i, i);
    miu_ssT.slice(vv) = miu_ssT_v;
    
    // accumulate for A M-step
    A_part1 += y_v.t() * miu_s_v.t(); // (nT x 1) * (1 x nL) -> nT x nL
    A_part2 += miu_ssT_v;             // nL x nL
    
    // keep last E_v_inv and Sigma_s_v (matching the R code which returns them from last vv)
    last_E_v_inv = E_v_inv;
    last_Sigma_s_v = Sigma_s_v;
  }
  
  // Build Estep list (matching R structure)
  List Estep = List::create(
    _["E_v_inv"] = wrap(last_E_v_inv),
    _["Sigma_s_v"] = wrap(last_Sigma_s_v),
    _["miu_s"] = wrap(miu_s),
    _["miu_ssT"] = wrap(miu_ssT),
    _["var_s"] = wrap(var_s)
  );
  
  if (return_MAP) {
    return Estep;
  }
  
  // Compute A_hat = A_part1 %*% solve(A_part2)
  arma::mat A_hat;
  // numerically safer to solve linear system for each column: A_hat = A_part1 * inv(A_part2)
  A_hat = A_part1 * inv(A_part2);
  
  // scale function: center columns, divide by column sds (matching R's scale())
  auto scale_columns = [&](arma::mat& M) {
    int nc = M.n_cols;
    for (int j = 0; j < nc; ++j) {
      arma::colvec col = M.col(j);
      double m = arma::mean(col);
      double s = arma::stddev(col, 1); // sample sd (default Armadillo uses N-1 denom)
      if (s == 0.0) s = 1.0;
      M.col(j) = (M.col(j) - m) / s;
    }
  };
  
  // If H provided, do Hinv %*% A_hat, scale, then H %*% A_hat
  if (H.isNotNull() && Hinv.isNotNull()) {
    arma::mat Hmat = as<arma::mat>(H);
    arma::mat Hinvmat = as<arma::mat>(Hinv);
    A_hat = Hinvmat * A_hat;
    scale_columns(A_hat);
    A_hat = Hmat * A_hat;
  } else {
    // no dimension reduction case: just scale
    scale_columns(A_hat);
  }
  
  // M-step for nu0_sq
  arma::mat Cinv = diagmat(1.0 / C_diag);        // nT x nT
  arma::mat Cinv_A = Cinv * A_hat;               // nT x nL
  arma::mat At_Cinv_A = A_hat.t() * Cinv * A_hat; // nL x nL
  
  double nu0sq_hat;
  if (update_nu0sq) {
    // nu0sq_part1 <- sum(diag(Cinv %*% t(BOLD) %*% BOLD))
    arma::mat BtB = BOLD.t() * BOLD;            // nT x nT
    double nu0sq_part1 = trace(Cinv * BtB);
    
    // nu0sq_part2 <- sum(diag(Cinv_A %*% t(miu_s) %*% BOLD))
    // t(miu_s) : nL x nV ; BOLD: nV x nT -> product: nL x nT ; Cinv_A (nT x nL) * (nL x nT) -> nT x nT
    arma::mat t_miu_s = miu_s.t();              // nL x nV
    arma::mat temp2 = Cinv_A * (t_miu_s * BOLD); // nT x nT
    double nu0sq_part2 = trace(temp2);
    
    // nu0sq_part3 <- sum(diag(At_Cinv_A %*% apply(miu_ssT, 2:3, sum) ))
    // compute sum over v of miu_ssT[v,,] -> nL x nL
    arma::mat sum_miu_ssT(nL, nL, fill::zeros);
    for (int vv = 0; vv < nV; ++vv) sum_miu_ssT += miu_ssT.slice(vv);
    double nu0sq_part3 = trace(At_Cinv_A * sum_miu_ssT);
    
    nu0sq_hat = (1.0 / ( (double)nT * (double)nV )) * (nu0sq_part1 - 2.0 * nu0sq_part2 + nu0sq_part3);
    
  } else {
    nu0sq_hat = as<double>(theta["nu0_sq"]);
  }
  
  // assemble theta_new list to match original R output
  List theta_new;
  theta_new["A"] = wrap(A_hat);
  theta_new["nu0_sq"] = wrap(nu0sq_hat);
  theta_new["Estep"] = Estep;
  
  // compute LL using R function compute_LL_std (must exist in R environment)
  // in R code: theta_new$LL[1] <- compute_LL_std(theta_new, prior_mean, prior_var, C_diag, BOLD, verbose=FALSE)
  try {
    Function compute_LL_std("compute_LL_std");
    NumericVector LLvec = compute_LL_std(theta_new, wrap(prior_mean), wrap(prior_var), wrap(C_diag), wrap(BOLD), Named("verbose") = false);
    // ensure LLvec has at least one element
    if (LLvec.size() > 0) {
      NumericVector LLslot = NumericVector::create(LLvec[0]);
      theta_new["LL"] = LLslot;
    } else {
      theta_new["LL"] = NumericVector::create(NA_REAL);
    }
  } catch(...) {
    // If compute_LL_std is not found or errors, silently set NA (to mimic safe behaviour)
    theta_new["LL"] = NumericVector::create(NA_REAL);
    if (verbose) Rcout << "Warning: compute_LL_std not found or error when calling it. theta_new$LL set to NA.\n";
  }
  
  return theta_new;
}