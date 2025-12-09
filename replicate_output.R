# required packages
#install.packages("Rcpp")
#install.packages("RcppArmadillo")
library(Rcpp)
library(RcppArmadillo)

# compile
# ---------------------------
# Dimensions for test example
# ---------------------------
nV <- 5   # number of voxels
nL <- 3   # number of latent components
nT <- 4   # number of time points


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

Rcpp::sourceCpp("/N/slate/nsagar/S610-FinalProject/updateTheta_BrainMap_independent.cpp")

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


