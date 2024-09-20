###################
### M-step ---- ###
###################

source("utils.R")

# param
# param:
param <- function(y, tau) {
  mu_new <- rep(0, ncol(tau))
  sigma_new <- rep(0, ncol(tau))
  for (g in 1:ncol(tau)) {
    mu_new[g] <- t(tau[,g]) %*% y / sum(tau[,g])
    sigma_new[g] <- t(tau[,g]) %*% (y - mu_new[g])^2 / sum(tau[,g])
  }
  return(list(mu=mu_new, sigma=sqrt(sigma_new)))
}

# @ gamma_f
gamma_fh <- function(X_tild, tau, eta_h, gamma_h, learning_rate=0.01) {
  for(g in 2:ncol(eta_h)){
    gradient <- colSums((tau[,g] - eta_h[,g])*X_tild)
    H <- matrix(0, nrow = ncol(X_tild), ncol = ncol(X_tild))
    for (i in 1:nrow(X_tild)) {
      H <- H - X_tild[i,]%*%t(X_tild[i,])*eta_h[i,g]*(1-eta_h[i,g])
    }
    H.inv <- solve(H)
    gamma_h[,g] <- gamma_h[,g] - learning_rate * H.inv %*% gradient
  }
  return(gamma_h)
}

gamma_f <- function(X_tild, tau, eta_h, gamma_h, learning_rate=0.005) {
  for(g in 2:ncol(eta_h)){
    gradient <- colSums((tau[,g] - eta_h[,g])*X_tild)
    gamma_h[,g] <- gamma_h[,g] + learning_rate * gradient
  }
  return(gamma_h)
}

# @
Mstep <- function(X, y, tau, eta, gamma_h) {
  param_new <- param(y, tau)
  mu_new <- param_new$mu
  sigma_new <- param_new$sigma
  gamma_new <- gamma_f(X, tau, eta, gamma_h)
  return(list(mu=mu_new, sigma=sigma_new, gamma=gamma_new))
}


# Standard model ----
param_stand <- function(X_tild, y, tau) {
  beta_new <- matrix(rep(0, ncol(X_tild)*ncol(tau)), ncol=ncol(tau))
  sigma_new <- rep(0, ncol(tau))
  for (g in 1:ncol(tau)) {
    T_g <- diag(tau[,g])
    beta_new[,g] <- solve(t(X_tild) %*% T_g %*% X_tild) %*% t(X_tild) %*% T_g %*% y
    sigma_new[g] <- t(tau[,g]) %*% (y - X_tild %*% beta_new[,g])^2 / sum(tau[,g])
  }
  return(list(beta=beta_new, sigma=sqrt(sigma_new)))
}

Mstep_stand <- function(X, y, tau, eta, gamma_h) {
  param_new <- param_stand(X, y, tau)
  beta_new <- param_new$beta
  sigma_new <- param_new$sigma
  gamma_new <- gamma_f(X, tau, eta, gamma_h)
  return(list(beta=beta_new, sigma=sigma_new, gamma=gamma_new))
}

