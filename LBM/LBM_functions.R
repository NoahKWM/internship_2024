source("Estep.R")
source("Mstep.R")
source("utils.R")

# LBM ----
lbm_f <- function(X, Q, L) {
  # init
  param <- init_param(X, Q, L)
  alpha_new <- param$alpha
  beta_new <- param$beta
  tau_new <- param$Tau_mat
  nu_new <- param$Nu_mat
  theta_new <- param$theta
  
  # E-step
  
  # M-step
}

lbm <- function(X, Q, L, max_iter=10) {
  param <- init_param(X, Q, L)
  alpha_new <- param$alpha
  beta_new <- param$beta
  tau_new <- param$Tau_mat
  nu_new <- param$Nu_mat
  theta_new <- param$theta
    
  for (iter in 1:max_iter) {
    # E-step
    E_step <- Estep(X, alpha_new, beta_new, theta_new, nu_new, tau_new)
    tau_new <- E_step$tau
    nu_new <- E_step$nu
    
    # M-step
    M_step <- Mstep(X, tau_new, nu_new)
    alpha_new <- M_step$alpha
    beta_new <- M_step$beta
    theta_new <- M_step$theta
    
  }
  
  return(list(
    alpha=alpha_new,
    beta=beta_new,
    theta=theta_new,
    tau=tau_new,
    nu=nu_new
  ))
}

lbm_bis <- function(X, Q, L, max_iter=100, tolerance=1e-6, nrep=5) {
  # Initialisation des tableaux pour stocker les résultats
  tab_log <- numeric(nrep)
  tab_alpha <- vector("list", nrep)
  tab_beta <- vector("list", nrep)
  tab_theta <- vector("list", nrep)
  tab_tau <- vector("list", nrep)
  tab_nu <- vector("list", nrep)
  
  for (k in 1:nrep) {
    # Init
    cat("Init numéro ", k, "\n")
    
    param <- init_param(X, Q, L)
    alpha_new <- param$alpha
    beta_new <- param$beta
    tau_new <- param$Tau_mat
    nu_new <- param$Nu_mat
    theta_new <- param$theta
    
    log_likelihood <- numeric(max_iter)
    
    for (iter in 1:max_iter) {
      # E-step
      E_step <- Estep(X, alpha_new, beta_new, theta_new, nu_new, tau_new)
      tau_new <- E_step$tau
      nu_new <- E_step$nu
      
      log_likelihood[iter] <- E_step$elbo
      
      if (any(is.na(log_likelihood[iter]))) {
        warning("NA value in log likelihood", iter)
        break
      }
      
      # M-step
      M_step <- Mstep(X, tau_new, nu_new)
      alpha_new <- M_step$alpha
      beta_new <- M_step$beta
      theta_new <- M_step$theta
      
      # Critère d'arrêt
      if (iter > 1) {
        if (abs(log_likelihood[iter] - log_likelihood[iter-1]) < tolerance) {
          print(paste("Convergence atteinte à l'itération", iter))
          break
        }
      }
    }
    # Stocker uniquement la dernière valeur de log_likelihood
    tab_log[k] <- log_likelihood[iter]
    tab_alpha[[k]] <- alpha_new
    tab_beta[[k]] <- beta_new
    tab_theta[[k]] <- theta_new
    tab_tau[[k]] <- tau_new
    tab_nu[[k]] <- nu_new
  }
  
  # Trouver l'index de la meilleure estimation basée sur log_likelihood
  best_index <- which.max(tab_log)
  
  # Retourner les meilleurs résultats
  return(invisible(list(
    best_log = tab_log[best_index],
    tau = tab_tau[[best_index]],
    nu = tab_nu[[best_index]],
    alpha = tab_alpha[[best_index]],
    beta = tab_beta[[best_index]],
    theta = tab_theta[[best_index]]
  )))
}


# co LBM -----
co_lbm <- function(X, Q, L, max_iter=100, tolerance=1e-6, nrep=5) {
  # Initialisation des tableaux pour stocker les résultats
  tab_log <- numeric(nrep)
  tab_alpha <- vector("list", nrep)
  tab_beta <- vector("list", nrep)
  tab_theta <- vector("list", nrep)
  tab_tau <- vector("list", nrep)
  tab_nu <- vector("list", nrep)
  
  for (k in 1:nrep) {
    # Init
    cat("Init numéro ", k, "\n")
    
    param <- co_init_param(X, Q, L)
    alpha_new <- param$alpha
    beta_new <- param$beta
    tau_new <- param$Tau_mat
    nu_new <- param$Nu_tensor
    theta_new <- param$theta
    
    log_likelihood <- numeric(max_iter)
    
    for (iter in 1:max_iter) {
      # E-step
      E_step <- co_Estep(X, alpha_new, beta_new, theta_new, nu_new, tau_new)
      tau_new <- E_step$tau
      nu_new <- E_step$nu
      
      log_likelihood[iter] <- E_step$elbo
      
      if (any(is.na(log_likelihood[iter]))) {
        warning("NA value in log likelihood", iter)
        break
      }
      
      # M-step
      M_step <- co_Mstep(X, tau_new, nu_new)
      alpha_new <- M_step$alpha
      beta_new <- M_step$betas
      theta_new <- M_step$theta
      
      # Critère d'arrêt
      if (iter > 1) {
        if (abs(log_likelihood[iter] - log_likelihood[iter-1]) < tolerance) {
          print(paste("Convergence atteinte à l'itération", iter))
          break
        }
      }
    }
    # Stocker uniquement la dernière valeur de log_likelihood
    tab_log[k] <- log_likelihood[iter]
    tab_alpha[[k]] <- alpha_new
    tab_beta[[k]] <- beta_new
    tab_theta[[k]] <- theta_new
    tab_tau[[k]] <- tau_new
    tab_nu[[k]] <- nu_new
  }
  
  # Trouver l'index de la meilleure estimation basée sur log_likelihood
  best_index <- which.max(tab_log)
  
  # Retourner les meilleurs résultats
  return(invisible(list(
    best_log = tab_log[best_index],
    tau = tab_tau[[best_index]],
    nu = tab_nu[[best_index]],
    alpha = tab_alpha[[best_index]],
    betas = tab_beta[[best_index]],
    theta = tab_theta[[best_index]]
  )))
}

