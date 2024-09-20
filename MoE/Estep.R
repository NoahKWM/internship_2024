###################
### E-step ---- ###
###################

source("utils.R")

# Simple model -----

# @ tau
# @ param: eta, mu (mean), sigma (écart-type), y (target var)
tau_f <- function(eta, mu, sigma, y) {
  res <- eta * sapply(1:ncol(eta), function(g) dnorm(y, mean = mu[g], sd = sigma[g]))
  res <- res/rowSums(res)
  return(res)
}

# @ ll
# @ param: tau, eta, mu, sigma, y
ll_f <- function(tau, eta, mu, sigma,y) {
  temp <- sum(tau * (log(.Machine$double.eps+eta) + sapply(1:ncol(eta), function(g) log(.Machine$double.eps+ dnorm(y, mu[g], sqrt(sigma[g]))))))
  return(temp)
}

# @ eta
# @ param: gamma, X
eta_f <- function(gamma, X) {
  lin_h <- X %*% gamma
  eta_h <- Softmax_log(lin_h)
  return(eta_h)
}

Estep <- function(X, y, gamma, mu, sigma) {
  eta_new <- eta_f(gamma, X)
  tau_new <- tau_f(eta_new, mu, sigma, y)
  ll <- ll_f(tau_new, eta_new, mu, sigma, y)
  return(list(tau = tau_new, ll = ll, eta = eta_new))
}


# Standard model ----

tau_f_stand <- function(X, y, beta_h, sigma, eta) {
  res <- eta * sapply(1:ncol(eta), function(g) dnorm(y, mean = X%*%beta_h[,g], sd = sigma[g]))
  res <- res/rowSums(res)
  return(res)
}

ll_f_stand <- function(X_tild, y, beta, sigma, tau, eta) {
  group <- apply(eta, 1, function(p) sample(1:ncol(tau), 1, prob = p))
  means <- numeric(length(y))
  for (i in 1:length(y)) {
    means[i] <- X_tild[i, ] %*% beta[, group[i]]
  }
  temp <- sum(tau * (log(.Machine$double.eps+eta) + sapply(1:ncol(eta), function(g) log(.Machine$double.eps+ dnorm(y, means, sqrt(sigma[g]))))))
  return(temp)
}

Estep_stand <- function(X, y, gamma, beta, sigma) {
  eta_new <- eta_f(gamma, X)
  tau_new <- tau_f_stand(X_tild, y, beta, sigma, eta_new)
  ll <- ll_f_stand(X_tild, y, beta, sigma, tau_new, eta_new)
  return(list(tau = tau_new, ll = ll, eta = eta_new))
}