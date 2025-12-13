# required packages
#install.packages("Rcpp")
#install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)

# compile
# ---------------------------
# Dimensions for test example
# ---------------------------
nV <- 50   # number of voxels
nL <- 30   # number of latent components
nT <- 40   # number of time points


prior_mean <- matrix(rnorm(nV * nL), nV, nL)

prior_var <- matrix(rexp(nV * nL, rate = 1), nV, nL)

BOLD <- matrix(rnorm(nV * nT), nV, nT)

theta <- list(
  A = matrix(rnorm(nT * nL), nT, nL),
  nu0_sq = 1.5   # any positive value
)

C_diag <- runif(nT, 0.5, 2)

H <- NULL
Hinv <- NULL

source("UpdateTheta_BrainMap.independent.R")
out <- UpdateTheta_BrainMap.independent(
  prior_mean = prior_mean,
  prior_var  = prior_var,
  BOLD       = BOLD,
  theta      = theta,
  C_diag     = C_diag,
  H          = H,
  Hinv       = Hinv,
  update_nu0sq = TRUE,
  return_MAP = FALSE,
  verbose = TRUE
)
out

Rcpp::sourceCpp("updateTheta_BrainMap_independent.cpp")

out2 <- UpdateTheta_BrainMap_independent_cpp(
  prior_mean = prior_mean,
  prior_var  = prior_var,
  BOLD       = BOLD,
  theta      = theta,
  C_diag     = C_diag,
  H          = H,
  Hinv       = Hinv,
  update_nu0sq = TRUE,
  return_MAP = FALSE,
  verbose = TRUE
)
out2

## =========================================================
## Validation tests: R vs C++ implementations
## =========================================================

tol <- 1e-8  # numerical tolerance

# Helper function for matrix comparison
compare_matrix <- function(M1, M2, name) {
  if (!all(dim(M1) == dim(M2))) {
    cat(name, ": DIMENSION MISMATCH\n")
    return(FALSE)
  }
  max_diff <- max(abs(M1 - M2))
  cat(sprintf("%s: max |diff| = %.3e\n", name, max_diff))
  max_diff < tol
}

# Helper function for scalar comparison
compare_scalar <- function(x1, x2, name) {
  diff <- abs(x1 - x2)
  cat(sprintf("%s: |diff| = %.3e\n", name, diff))
  diff < tol
}

# Helper function for 3D array comparison
compare_array3 <- function(A1, A2, name) {
  if (!all(dim(A1) == dim(A2))) {
    cat(name, ": DIMENSION MISMATCH\n")
    return(FALSE)
  }
  max_diff <- max(abs(A1 - A2))
  cat(sprintf("%s: max |diff| = %.3e\n", name, max_diff))
  max_diff < tol
}

## -------------------------
## 1. Compare A
## -------------------------
A_ok <- compare_matrix(out$A, out2$A, "A")

## -------------------------
## 2. Compare nu0_sq
## -------------------------
nu0_ok <- compare_scalar(out$nu0_sq, out2$nu0_sq, "nu0_sq")

## -------------------------
## 3. Compare Estep components
## -------------------------

E_v_inv_ok <- compare_matrix(
  out$Estep$E_v_inv,
  out2$Estep$E_v_inv,
  "Estep$E_v_inv"
)

Sigma_ok <- compare_matrix(
  out$Estep$Sigma_s_v,
  out2$Estep$Sigma_s_v,
  "Estep$Sigma_s_v"
)

miu_s_ok <- compare_matrix(
  out$Estep$miu_s,
  out2$Estep$miu_s,
  "Estep$miu_s"
)

miu_ssT_ok <- compare_array3(
  out$Estep$miu_ssT,
  out2$Estep$miu_ssT,
  "Estep$miu_ssT"
)

var_s_ok <- compare_matrix(
  out$Estep$var_s,
  out2$Estep$var_s,
  "Estep$var_s"
)

## -------------------------
## 4. Final summary
## -------------------------
cat("\n====== SUMMARY ======\n")
results <- c(
  A = A_ok,
  nu0_sq = nu0_ok,
  E_v_inv = E_v_inv_ok,
  Sigma_s_v = Sigma_ok,
  miu_s = miu_s_ok,
  miu_ssT = miu_ssT_ok,
  var_s = var_s_ok
)

print(results)

if (all(results)) {
  cat("SUCCESS: All outputs match within tolerance.\n")
} else {
  cat("WARNING: Some outputs do NOT match.\n")
}

