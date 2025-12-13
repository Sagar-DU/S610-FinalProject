UpdateTheta_BrainMap.independent <- function(
    prior_mean, prior_var, BOLD,
    theta, C_diag, H, Hinv,
    update_nu0sq =TRUE, return_MAP=FALSE, verbose=TRUE){

  nL <- ncol(prior_mean)
  nV <- nrow(BOLD)
  nT <- ncol(BOLD)

  #initialize new objects
  theta_new <- list(A = matrix(NA, nT, nL), nu0_sq = NA)
  #two parts of product for A-hat (construct each looping over voxels)
  A_part1 <- matrix(0, nT, nL)
  A_part2 <- matrix(0, nL, nL)

  A <- theta$A # TxQ
  nu0_sq <- theta$nu0_sq # const
  nu0C_inv <- diag(1/(C_diag*nu0_sq)) #Sigma0_inv in matlab code # TxT
  At_nu0Cinv <- crossprod(A, nu0C_inv) # QxT
  At_nu0Cinv_A <- At_nu0Cinv %*% A # QxQ

  #store posterior moments for M-step of nu0_sq
  miu_s <- matrix(NA, nrow=nV, ncol=nL)
  miu_ssT <- array(NA, dim=c(nV, nL, nL))
  var_s <- matrix(NA, nrow=nV, ncol=nL)

  for (vv in 1:nV) {
    y_v <- BOLD[vv,] # T
    s0_v <- prior_mean[vv,] # Q

    ##########################################
    ### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
    ##########################################

    E_v_inv <- diag(1/prior_var[vv,]) # QxQ
    Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A) # QxQ
    miu_s_v <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    miu_ssT_v <- tcrossprod(miu_s_v) + Sigma_s_v #QxQ
    miu_s[vv,] <- miu_s_v #save for M-step of nu0_sq
    miu_ssT[vv,,] <- miu_ssT_v #save for M-step of nu0_sq
    var_s[vv,] <- diag(Sigma_s_v)

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 <- A_part1 + tcrossprod(as.matrix(y_v), miu_s_v) #TxQ + (Tx1)*(1xQ)
    A_part2 <- A_part2 + miu_ssT_v #QxQ
  }

  Estep <- list(
    E_v_inv=E_v_inv,
    Sigma_s_v=Sigma_s_v,
    miu_s=miu_s,
    miu_ssT=miu_ssT,
    var_s=var_s
  )

  if(return_MAP) return(Estep)

  #A_hat <- orthonorm(A_part1 %*% solve(A_part2))
  A_hat <- (A_part1 %*% solve(A_part2))
  #fix scale of A
  if(!is.null(H)){ #case where dimension reduction is in place
    A_hat <- Hinv %*% A_hat
    A_hat <- scale(A_hat)
    A_hat <- H %*% A_hat
  } else { #case with no dimension reduction
    A_hat <- scale(A_hat)
  }

  ##########################################
  ### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
  ##########################################

  #cat('Updating Error Variance nu0_sq \n')

  Cinv <- diag(1/C_diag)
  Cinv_A <- Cinv %*% A_hat
  At_Cinv_A <- t(A_hat) %*% Cinv %*% A_hat

  if (update_nu0sq) {
    # #old version
    # nu0sq_part1 <- nu0sq_part2 <- nu0sq_part3 <- 0
    # for(vv in 1:nV){
    #
    #   y_v <- BOLD[vv,]
    #   nu0sq_part1 <- nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
    #   nu0sq_part2 <- nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[vv,]
    #   nu0sq_part3 <- nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[vv,,]))
    # }
    #new version (use Trace instead of summing over v, and use trick that Tr(AB) = TR(BA) when AB is VxV but BA is QxQ or TxT)
    nu0sq_part1 <- sum(diag(Cinv %*% t(BOLD) %*% BOLD))
    nu0sq_part2 <- sum(diag(Cinv_A %*% t(miu_s) %*% BOLD))
    nu0sq_part3 <- sum(diag(At_Cinv_A %*% apply(miu_ssT, 2:3, sum) ))
    nu0sq_hat <- 1/(nT*nV)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)
  } else {
    #keep initial estimate of nu0sq
    nu0sq_hat <- theta$nu0_sq
  }

  # RETURN NEW PARAMETER ESTIMATES
  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  theta_new$Estep <- Estep

  # COMPUTE LL FOR SQUAREM
#  theta_new$LL[1] <- compute_LL_std(
 #   theta_new, prior_mean, prior_var, C_diag, BOLD, verbose=FALSE
  #)

  return(theta_new)
}