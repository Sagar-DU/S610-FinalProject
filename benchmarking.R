# install.packages("microbenchmark")
library(microbenchmark)
library(Rcpp)

source("UpdateTheta_BrainMap.independent.R")
Rcpp::sourceCpp("updateTheta_BrainMap_independent.cpp")

microbenchmark(
  UpdateTheta_BrainMap.independent(
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
  ),
  
  UpdateTheta_BrainMap_independent_cpp(
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
)