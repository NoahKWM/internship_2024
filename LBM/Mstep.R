source('utils.R')

# alpha ----
alpha_f <- function(tau) {
  vect <- rep(1, dim(tau)[1])
  alpha <- vect %*% tau / dim(tau)[1]
  return(alpha)
}

# beta ----
beta_f <- function(nu) {
  vect <- rep(1, dim(nu)[1])
  beta <- vect %*% nu / dim(nu)[1]
  return(beta)
}

# Theta ----
theta_f <- function(X, tau, nu) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ones <- matrix(1, nrow = n, ncol = p)
  den_mat <- t(tau)%*% ones %*% nu
  num_mu_mat <- t(tau)%*% X %*% nu
  
  Q <- dim(tau)[2]
  L <- dim(nu)[2]
  res <- array(NA, c(Q, L, 2))
  
  for (k in 1:Q) {
    for (l in 1:L) {
      res[k,l,1] <- num_mu_mat[k,l] / den_mat[k,l]
      res[k,l,2] <- sqrt((t(tau)[k,]%*%(X - res[k,l,1])^2 %*% nu[,l])/den_mat[k,l])
    }
  }
  return(res)
}

# M-step ----
Mstep <- function(X, tau, nu) {
  alpha <- alpha_f(tau)
  beta <- beta_f(nu)
  theta <- theta_f(X, tau, nu)
  return(list(alpha=alpha, beta=beta, theta=theta))
}

# ---------------

# alpha ----
co_alpha_f <- function(tau) {
  vect <- rep(1, dim(tau)[1])
  alpha <- vect %*% tau / dim(tau)[1]
  return(alpha)
}

# beta ----
co_beta_f <- function(nu) {
  betas <- matrix(NA, dim(nu)[3], dim(nu)[2])
  for(k in 1:dim(nu)[3]) {
    vect <- rep(1, dim(nu)[1])
    betas[k,] <- vect %*% nu[,,k] / dim(nu)[1]
  }
  return(betas)
}

# Theta ----
co_theta_f <- function(X, tau, nu) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Q <- dim(nu)[3]
  L <- dim(nu)[2]
  
  ones <- matrix(1, nrow = n, ncol = p)
  res <- array(NA, c(Q, L, 2))
  
  for(k in 1:Q) {
    den_mat <- t(tau)%*% ones %*% nu[,,k]
    num_mu_mat <- t(tau)%*% X %*% nu[,,k]
    
    for(l in 1:L) {
      res[k,l,1] <- num_mu_mat[k,l] / den_mat[k,l]
      res[k,l,2] <- sqrt((t(tau)[k,]%*%(X - res[k,l,1])^2 %*% nu[,l,k])/den_mat[k,l])
    }
  }
  
  return(res)
}

# M-step ----
co_Mstep <- function(X, tau, nu) {
  alpha <- co_alpha_f(tau)
  betas <- co_beta_f(nu)
  theta <- co_theta_f(X, tau, nu)
  return(list(alpha=alpha, betas=betas, theta=theta))
}