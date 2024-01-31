library(Matrix)
library(kernlab)
library(osqp)

rm(list = ls())
source("~/Documents/kernel-maic.R")

gen_data <- function(n, sig0 = 2, sig1 = 1.5,
                     scenario = c("anchor", "unanchor",
                                  "ps-mis-anchor", "out-mis-anchor",
                                  "ps-mis-unanchor", "out-mis-unanchor",)){
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp((x1 + x4)/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1))))
  u3 <- as.numeric(scale(log(abs(x2*x3))))
  u4 <- as.numeric(scale((x3 + x4)^2))
  
  # v1 <- as.numeric(scale(exp(x1/2)))
  # v2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  # v3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  # v4 <- as.numeric(scale((x2 + x4 + 10)^2))
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  # V <- cbind(int = rep(1, n), v1, v2, v3, v4)
  
  # coefficients
  beta0 <- c(2, -3, -1, 1, 3)
  beta1 <- c(0, 2, -2, -2, 2)
  alpha0 <- c(-2, -1, 3, -3, 1)
  alpha1 <- c(1, 2, -2, 2, -2)
  delta0 <- c(-0.25, 0, 0, -0, -0)
  delta1 <- c(0.25, 0, -0, 0, -0)
  gamma <- c(0, -0.5, 0.5, -0.5, 0.5)
  
  if (scenario == "anchor"){
    
    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-X %*% delta1)) + (1 - s)/(1 + exp(-X %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(X %*% beta0)
    mu_1 <- mu_0 + c(s*(X %*% alpha1) + (1 - s)*(X %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(X[s == 0,] %*% (alpha0 - alpha1))
    
  } else if (scenario == "unanchor"){
    
    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-X %*% delta1)) + (1 - s)/(1 + exp(-X %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(s*(X %*% beta1) + (1 - s)*(X %*% beta0))
    mu_1 <- mu_0 + c(s*(X %*% alpha1) + (1 - s)*(X %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(X[s == 0,] %*% (alpha0 - alpha1))
    
  } else if (scenario == "ps-mis-anchor") {
    
    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-U %*% delta1)) + (1 - s)/(1 + exp(-U %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(X %*% beta0)
    mu_1 <- mu_0 + c(s*(X %*% alpha1) + (1 - s)*(X %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(X[s == 0,] %*% (alpha0 - alpha1))
    
  } else if (scenario == "out-mis-anchor"){

    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-X %*% delta1)) + (1 - s)/(1 + exp(-X %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + c(s*(U %*% alpha1) + (1 - s)*(U %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(U[s == 0,] %*% (alpha0 - alpha1))
    
  } else if (scenario == "ps-mis-unanchor") {
    
    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-U %*% delta1)) + (1 - s)/(1 + exp(-U %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(s*(X %*% beta1) + (1 - s)*(X %*% beta0))
    mu_1 <- mu_0 + c(s*(X %*% alpha1) + (1 - s)*(X %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(X[s == 0,] %*% (alpha0 - alpha1))
    
  } else if (scenario == "out-mis-unanchor"){
    
    # sample score
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    
    # propensity score
    e_X <- c(s/(1 + exp(-X %*% delta1)) + (1 - s)/(1 + exp(-X %*% delta0)))
    z <- rbinom(n, 1, e_X) # treatment
    
    # outcome models
    mu_0 <- c(s*(U %*% beta1) + (1 - s)*(U %*% beta0))
    mu_1 <- mu_0 + c(s*(U %*% alpha1) + (1 - s)*(U %*% alpha0))
    
    # indirect treatment effect
    tau <- mean(U[s == 0,] %*% (alpha0 - alpha1))
    
  }
  
  # potential outcomes
  y_tmp0 <- rnorm(n, mu_0, sig0)
  y_tmp1 <- rnorm(n, mu_1, sig1)
  y_pot <- cbind(y_tmp0, y_tmp1)
  
  # observed outcome
  y <- z*y_pot[,2] + (1 - z)*y_pot[,1]
  
  # create simulation dataset
  sim <- list(y = y, z = z, s = s, X = X, U = U, tau = tau)
  
  return(sim)
  
}

# Fits the RKHS balancing weights using a variety of kernels
sim_fit <- function(idx = 1, simDat, ...) {

  dat <- simDat[,idx]
  tau <- dat$tau
  
  # data components
  Y <- dat$y
  Z <- dat$z
  S <- dat$s
  X <- dat$X
  
  # stratify by study
  Y0 <- Y[S == 0]
  Z0 <- Z[S == 0]
  X0 <- X[S == 0,]
  Y1 <- Y[S == 1]
  Z1 <- Z[S == 1]
  X1 <- X[S == 1,]
  
  # dimensions
  m <- nrow(X0)
  n <- nrow(X1)
  
  # linear kernel
  linear0 <- qpmaic(Y = Y0, X = X0, Z = Z0, target = X0, 
                    kernel = kernlab::vanilladot(), lambda = 0, 
                    eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE)
  linear1 <- qpmaic(Y = Y1, X = X1, Z = Z1, target = X0,
                    kernel = kernlab::vanilladot(), lambda = 0,
                    eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE)
  
  tau_linear <- linear0$mu1 - linear1$mu1
  # pool_linear <- c((m - 1)*linear0$sig1 + (n - 1)*linear1$sig1)/(n + m - 2)
  omega_linear <- linear0$sig1 + linear1$sig1
  
  # RBF kernel
  rbf0 <- qpmaic(Y = Y0, X = X0, Z = Z0, target = X0,
                 kernel = kernlab::rbfdot(), lambda = 1,
                 eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE)
  rbf1 <- qpmaic(Y = Y1, X = X1, Z = Z1, target = X0,
                 kernel = kernlab::rbfdot(), lambda = 1,
                 eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE)
  
  tau_rbf <- rbf0$mu1 - rbf1$mu1
  # pool_rbf <- c((m - 1)*rbf0$sig1 + (n - 1)*rbf1$sig1)/(n + m - 2)
  omega_rbf <- rbf0$sig1 + rbf1$sig1
  
  # combine results
  est <- c(tau = tau, linear = tau_linear, rbf = tau_rbf)
  se <- c(linear = sqrt(omega_linear), rbf = sqrt(omega_rbf))
  
  return(list(est = est, se = se))
  
}

n.iter <- 100 # simulation iterations

# run simulation study
simDat <- replicate(100, gen_data(n = 2000, scenario = "out-mis-anchor"))
simFit <- lapply(1:100, sim_fit, simDat = simDat)

# summarize
est_mat <- do.call(rbind, lapply(simFit, function(lst, ...) lst$est))
se_mat <- do.call(rbind, lapply(simFit, function(lst, ...) lst$se))
cover_mat <- cbind(linear = est_mat[,2] - 1.96*se_mat[,1] <= mean(est_mat[,1]) &
                     est_mat[,2] + 1.96*se_mat[,1] >= mean(est_mat[,1]),
                   rbf = est_mat[,3] - 1.96*se_mat[,2] <= mean(est_mat[,1]) &
                     est_mat[,3] + 1.96*se_mat[,2] >= mean(est_mat[,1]))

# results
colMeans(abs(est_mat[,2:3] - mean(est_mat[,1])))
colMeans(se_mat)
colMeans(cover_mat)
