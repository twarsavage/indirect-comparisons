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
qpmaic <- function(Y, X, Z, target, kernel = kernlab::vanilladot(), lambda = 0,
                   eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE, ...) {
  
  # ensure that covariate matrices are matrices and get total number of units
  X <- as.matrix(X)
  m <- nrow(target)
  n <- nrow(X)
  
  # compute number of treated and control units within each site
  n1 <- sum(Z)
  n0 <- sum(1 - Z)
  
  # construct linear term vector
  Q <- kernlab::kernelMatrix(kernel = kernel, x = X, y = target)
  q <- -c(m*rowMeans(Q)) # maybe use n and pooled se?
  
  # Compute kernel matrices (Z = 1)
  kern1 <- kernlab::kernelMatrix(kernel, X)
  P1 <- tcrossprod(Z) * kern1
  
  # Compute kernel matrices (Z = 0)
  kern0 <- kernlab::kernelMatrix(kernel, X)
  P0 <- tcrossprod(1 - Z) * kern0
  
  # construct quadratic matrix
  P <- P0 + P1 + lambda * Matrix::Diagonal(n)
  # P <- m*(P0/n0 + P1/n1 + lambda * Matrix::Diagonal(n, c(Z/n1 + (1 - Z)/n0)))/2
  

  # sum to n weights
  A1 <- Matrix::rbind2(Matrix::t(Z), Matrix::t(1 - Z))
  l1 <- c(m, m) # maybe use n and pooled se?
  u1 <- c(m, m)
  
  # upper and lower bounds of individual weights
  A2 <- Matrix::Diagonal(n)
  l2 <- rep(0, n)
  u2 <- rep(m, n) # maybe use n and pooled se?
  
  A <- rbind(A1, A2)
  l <- c(l1, l2)
  u <- c(u1, u2)
  
  # set optimization settings
  settings <- do.call(osqp::osqpSettings,
                      c(list(verbose = verbose,
                             eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))
  
  # solve optimization problem
  solution <- osqp::solve_osqp(P = P, q = q, A = A, l = l, u = u, pars = settings)
  weights <- solution$x
  
  # compute imbalances
  imbalance1 <- c(colMeans(target)) - c(t(Z*X) %*% weights)/m
  names(imbalance1) <- colnames(X)
  imbalance0 <- c(colMeans(target)) - c(t((1 - Z)*X) %*% weights)/m
  names(imbalance0) <- colnames(X)
  
  # means and effect estimates
  mu1 <- sum(weights*Z*Y)/m
  mu0 <- sum(weights*(1 - Z)*Y)/m
  tau <- mu1 - mu0
  
  # variances (needs work)
  sig1 <- var(weights*Z*(Y - mu1))/m
  sig0 <- var(weights*(1 - Z)*(Y - mu0))/m
  omega <- sig1 + sig0
  
  # return output
  return(list(tau = tau, omega = omega, mu1 = mu1, mu0 = mu0, sig1 = sig1, sig0 = sig0,
              weights = weights, imbalance0 = imbalance0, imbalance1 = imbalance1))
  
}

