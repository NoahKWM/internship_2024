source("utils.R")

# tau ----
log_Tik_f <- function(X, alpha, theta, nu, i, k) {
  log_T_ik <- -1
  for (j in 1:dim(nu)[1]) {
    for (l in 1:dim(nu)[2]) {
      log_T_ik <- log_T_ik + nu[j,l]*log(.Machine$double.eps + dnorm(X[i,j], mean = theta[k,l,1], sd = .Machine$double.eps +theta[k,l,2])) + log(alpha[k])
    }
  }
  return(log_T_ik)
}

tau_f <- function(X, alpha, theta, nu, Q) {
  tmp <- matrix(0, nrow = dim(X)[1], ncol = Q)
  for (i in 1:dim(X)[1]) {
    for (k in 1:Q) {
      tmp[i,k] <- log_Tik_f(X, alpha, theta, nu, i, k)
    }
  }
  tau <- Softmax_log(tmp)
  return(tau)
}

# nu ----
log_Vjl_f <- function(X, beta, theta, tau, j, l){
  log_V_jl <- -1
  for (i in 1:dim(tau)[1]) {
    for (k in 1:dim(tau)[2]) {
      log_V_jl <- log_V_jl + log(.Machine$double.eps+dnorm(X[i,j], mean = theta[k,l,1], sd = .Machine$double.eps + theta[k,l,2]))*tau[i,k] + log(beta[l])
    }
  }
  return(log_V_jl)
}

nu_f <- function(X, beta, theta, tau, L) {
  tmp <- matrix(0, nrow = dim(X)[2], ncol = L)
  for (j in 1:dim(X)[2]) {
    for (l in 1:L) {
      tmp[j,l] <- log_Vjl_f(X, beta, theta, tau, j, l)
    }
  }
  nu <- Softmax_log(tmp)
  return(nu)
}

# elbo ----
elbo_f <- function(X, tau, nu, alpha, beta,theta) {
  ll <- 0
  for(i in 1:dim(X)[1]) {
    for(j in 1:dim(X)[2]) {
      for(k in 1:dim(tau)[2]) {
        ll <- ll + tau[i,k]*log(alpha[k]) - tau[i,k]*log(tau[i,k])
        for(l in 1:dim(nu)[2]) {
          ll <- ll + tau[i,k]*nu[j,l]*log(dnorm(X[i,j], mean = theta[k,l,1], sd = .Machine$double.eps + theta[k,l,2]))
          ll <- ll + nu[j,l]*log(beta[k]) - nu[j,l]*log(nu[j,l])
        }
      }
    }
  }
  return(ll)
}

# E-step ----
Estep <- function(X, alpha, beta, theta, nu, tau) {
  Q <- dim(theta)[1]
  L <- dim(theta)[2]
  tau_new <- tau_f(X, alpha, theta, nu, Q)
  nu_new <- nu_f(X, beta, theta, tau_new, L)
  elbo <- elbo_f(X, tau, nu, alpha, beta, theta) 
  return(list(tau=tau_new, nu=nu_new,elbo=elbo))
  #return(list(tau=tau_new, nu=nu_new))
}

# -------------------------


# tau ----
co_log_Tik_f <- function(X, alpha, theta, nu, i, k) {
  log_T_ik <- -1
  for (j in 1:(dim(nu)[1])) {
    for (l in 1:(dim(nu)[2])) {
      log_T_ik <- log_T_ik + nu[j,l,k]*log(.Machine$double.eps + dnorm(X[i,j], mean = theta[k,l,1],
                                                                       sd = .Machine$double.eps + theta[k,l,2])) + log(alpha[k])
    }
  }
  return(log_T_ik)
}

co_tau_f <- function(X, alpha, theta, nu, Q) {
  tmp <- matrix(NA, nrow = dim(X)[1], ncol = Q)
  for (i in 1:(dim(X)[1])) {
    for (k in 1:Q) {
      tmp[i,k] <- co_log_Tik_f(X, alpha, theta, nu, i, k)
    }
  }
  tau <- Softmax_log(tmp)
  return(tau)
}

# nu ----
co_log_Vkjl_f <- function(X, beta, theta, tau, k, j, l){
  log_V_jl <- -1
  for (i in 1:dim(tau)[1]) {
    log_V_jl <- log_V_jl + log(.Machine$double.eps+dnorm(X[i,j], mean = theta[k,l,1], sd = theta[k,l,2]))*tau[i,k] + log(beta[l])
  }
  return(log_V_jl)
}

co_nu_f <- function(X, beta, theta, tau, L) {
  nu <- array(0, dim=c(dim(X)[2], L, dim(tau)[2]))
  for(k in 1:dim(tau)[2]) {
    tmp <- matrix(0, nrow = dim(X)[2], ncol = L)
    for (j in 1:dim(X)[2]) {
      for (l in 1:L) {
        tmp[j,l] <- co_log_Vkjl_f(X, beta, theta, tau, k, j, l)
      }
    }
    nu[,,k] <- Softmax_log(tmp)
  }
  return(nu)
}

# elbo
co_elbo_f <- function(X, tau, nu, alpha, beta,theta) {
  ll <- 0
  for(i in 1:dim(X)[1]) {
    for(j in 1:dim(X)[2]) {
      for(k in 1:dim(tau)[2]) {
        ll <- ll + tau[i,k]*log(alpha[k]) - tau[i,k]*log(tau[i,k])
        for(l in 1:dim(nu)[2]) {
          ll <- ll + tau[i,k]*nu[j,l,k]*log(dnorm(X[i,j], mean = theta[k,l,1], sd = .Machine$double.eps + theta[k,l,2]))
          ll <- ll + nu[j,l,k]*log(beta[k,l]) - nu[j,l,k]*log(nu[j,l,k])
        }
      }
    }
  }
  return(ll)
}

# E-step
co_Estep <- function(X, alpha, beta, theta, nu, tau) {
  Q <- dim(theta)[1]
  L <- dim(theta)[2]
  tau_new <- co_tau_f(X, alpha, theta, nu, Q)
  nu_new <- co_nu_f(X, beta, theta, tau_new, L)
  elbo <- co_elbo_f(X, tau, nu, alpha, beta, theta) 
  return(list(tau=tau_new, nu=nu_new,elbo=elbo))
  #return(list(tau=tau_new, nu=nu_new))
}