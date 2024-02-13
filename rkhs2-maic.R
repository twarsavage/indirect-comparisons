# Quadratic Program for matching adjusted indirect comparisons

#' Quadratic program for nonparametric matching adjusted indirect comparisons
#' @param Y vector of outcomes
#' @param Z vector of treatment assignment
#' @param X n x d matrix of covariates
#' @param target n x d matrix of target covariate distribution
#' @param lambda Regularization hyper parameter, default 0
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#' @export
qp_maic <- function(Y, X, Z, S, kernel = kernlab::vanilladot(), lambda = 0,
                    pscores = 1, eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE, ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  X <- as.matrix(X)
  m <- sum(1 - S)
  n <- sum(S)
  
  # idx <- apply(X, 2, function(vec) {
  #   if (is.numeric(vec)){
  #     return(var(vec) != 0)
  #   } else
  #     return(FALSE)
  # })
  # 
  # X <- X[,which(idx)]
  
  # stratify by study
  Y0 <- as.numeric(Y[S == 0])
  Z0 <- as.numeric(Z[S == 0])
  X0 <- as.matrix(X[S == 0,])
  Y1 <- as.numeric(Y[S == 1])
  Z1 <- as.numeric(Z[S == 1])
  X1 <- as.matrix(X[S == 1,])
  
  # Propensity Score Estimates
  pscores0 <- c(SuperLearner(Y = Z0, X = data.frame(X0[,-1]), 
                             SL.library = c("SL.mean", "SL.glm", "SL.earth"),
                             family = binomial())$SL.predict)
  pscores1 <- c(SuperLearner(Y = Z1, X = data.frame(X1[,-1]), 
                             SL.library = c("SL.mean", "SL.glm", "SL.earth"),
                             family = binomial())$SL.predict)
  
  pscores <- rep(1, nrow(X))
  pscores[S == 0] <- pscores0
  pscores[S == 1] <- pscores1
  
  fitw <- standardize_treatment_kernel(X0 = X, Xtau = X, Xtarget = X0, S = S, Z = Z, pscores = pscores,
                                       kernel0 = kernel, kerneltau = kernel, lambda = lambda,
                                       lowlim = 0, uplim = m, scale_sample_size = FALSE,
                                       verbose = verbose, return_program = FALSE, init_uniform = FALSE,
                                       eps_abs = eps_abs, eps_rel = eps_rel, gc = FALSE)
  
  weights <- fitw$weights/(Z*pscores + (1 - Z)*(1 - pscores))
  weights0 <- weights[S == 0]
  weights1 <- weights[S == 1]
  
  # compute imbalances
  imbalance11 <- c(colMeans(X0)) - c(t(Z1*X1) %*% weights1)/n
  names(imbalance11) <- colnames(X)
  imbalance10 <- c(colMeans(X0)) - c(t((1 - Z1)*X1) %*% weights1)/n
  names(imbalance10) <- colnames(X)
  imbalance01 <- c(colMeans(X0)) - c(t(Z0*X0) %*% weights0)/m
  names(imbalance01) <- colnames(X)
  imbalance00 <- c(colMeans(X0)) - c(t((1 - Z0)*X0) %*% weights0)/m
  names(imbalance00) <- colnames(X)
  
  # means and effect estimates (S == 1)
  mu11 <- sum(weights1*Z1*Y1)/sum(Z1*weights1)
  eif11 <- weights1*Z1*Y1 - mu11
  sig11 <- var(eif11)/n
  
  mu10 <- sum(weights1*(1 - Z1)*Y1)/sum((1 - Z1)*weights1)
  eif10 <- weights1*(1 - Z1)*Y1 - mu10
  sig10 <- var(eif10)/n
  
  tau1 <- mu11 - mu10
  eif1 <- weights1*Z1*Y1 - weights1*(1 - Z1)*Y1 - tau1
  omega1 <- var(eif1)/n
  
  # means and effect estimates S == 0
  mu01 <- sum(weights0*Z0*Y0)/sum(Z0*weights0)
  eif01 <- weights0*Z0*Y0 - mu01
  sig01 <- var(eif01)/m
  
  mu00 <- sum(weights0*(1 - Z0)*Y0)/sum((1 - Z0)*weights0)
  eif00 <- weights0*(1 - Z0)*Y0 - mu00
  sig00 <- var(eif01)/m
  
  tau0 <- mu01 - mu00
  eif0 <- weights0*Z0*Y0 - weights0*(1 - Z0)*Y0 - tau0
  omega0 <- var(eif0)/m
  
  tau <- mu11 - mu01
  omega <- c((n - 1)*sig11 + (m - 1)*sig01) / (n + m - 2)
  
  # return output
  return(list(tau = tau, omega = omega, 
              tau0 = tau0, omega0 = omega0, 
              tau1 = tau1, omega1 = omega1,
              imbalance11 = imbalance11, imbalance10 = imbalance10,
              imbalance01 = imbalance01, imbalance00 = imbalance00))
  
}

# Quadratic Program for matching adjusted indirect comparisons

#' Quadratic program for nonparametric matching adjusted indirect comparisons
#' @param Y vector of outcomes
#' @param Z vector of treatment assignment
#' @param X n x d matrix of covariates
#' @param target n x d matrix of target covariate distribution
#' @param lambda Regularization hyper parameter, default 0
#' @param eps_abs Absolute error tolerance for solver
#' @param eps_rel Relative error tolerance for solver
#' @param ... Extra arguments for osqp solver
#' @export
drqp_maic <- function(Y, X, Z, target, kernel = kernlab::vanilladot(), lambda = 0,
                      family = gaussian(), sl.lib = c("SL.mean", "SL.glm", "SL.earth"), 
                      eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE, ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  X <- as.matrix(X)
  m <- nrow(target)
  n <- nrow(X)
  
  idx <- apply(X, 2, function(vec) {
    if (is.numeric(vec)){
      return(var(vec) != 0)
    } else
      return(FALSE)
  })
  
  X <- X[,which(idx)] # remove intercepts
  
  qp <- qp_maic(Y = Y, X = X, Z = Z, target = target, kernel = kernel, lambda = lambda,
                eps_abs = eps_abs, eps_rel = eps_rel, verbose = verbose) 
  
  weights <- qp$weights
  
  # simple outcome model
  mumod <- SuperLearner(Y = Y, X = data.frame(X, Z = Z), family = family, SL.library = sl.lib)
  muhat <- mumod$SL.predict
  muhat0 <- suppressWarnings(c(predict(mumod, newdata = data.frame(target, Z = 0), type = "response")$pred))
  muhat1 <- suppressWarnings(c(predict(mumod, newdata = data.frame(target, Z = 1), type = "response")$pred))
  
  # Point Estimates
  psi1 <- Z*c(Y - muhat)*weights
  psi0 <- (1 - Z)*c(Y - muhat)*weights 
  mu1 <- (sum(psi1) + sum(muhat1))/m
  mu0 <- (sum(psi0) + sum(muhat0))/m
  tau_eif <- 
    
    
    # variances (needs work)
    mu1_eif <- c(psi1, muhat1) - mu1
  mu0_eif <- c(psi0, muhat0) - mu0
  tau <- mu1 - mu0
  sig1 <- var(mu1_eif)/m
  sig0 <- var(mu0_eif)/m
  omega <- var(c(psi1, muhat1) - c(psi0, muhat0))/m
  
  # return output
  return(list(tau = tau, omega = omega, mu1 = mu1, mu0 = mu0, sig1 = sig1, sig0 = sig0,
              mu1_eif = mu1_eif, mu0_eif = mu0_eif, weights = weights, 
              imbalance0 = qp$imbalance0, imbalance1 = qp$imbalance1))
  
}
