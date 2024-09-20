Softmax_log <- function(log_X){
  if(!is.matrix(log_X)){log_X <- as.matrix(log_X)}
  K <- ncol(log_X)
  
  log_X <- log_X - apply(log_X,1,max)
  
  ## Now going back to exponential with the same normalization
  X <- exp(log_X) #(matrix(1,n,1) %*% pi) * exp(logX)
  X <- pmin(X,.Machine$double.xmax)
  X <- pmax(X,.Machine$double.xmin)
  X <- X / (rowSums(X) %*% matrix(1,1,K))
  X <- pmin(X,1-.Machine$double.eps)
  X <- pmax(X,.Machine$double.eps)
  
  return(X)
}

init_func <- function(G) {
  list(
    gamma = cbind(matrix(0, 3, 1), matrix(runif(3, -1, 1), 3, G-1)),
    mu = runif(G, -5, 5),                                           
    sigma = runif(G, 0, 2)                                          
  )
}

init_param <- function(X_tild, y, G) {
  res <- kmeans(X_tild[,-1], centers = G, nstart = 25)
  
  # Moyenne des groupes 
  mu_hat <- res$centers
  
  # Spooled
  cov_matrices <- list()
  group_sizes <- res$size
  for (cluster in unique(res$cluster)) {
    cluster_data <- X_tild[res$cluster == cluster, -1]
    cov_matrix <- cov(cluster_data)
    cov_matrices[[paste("Cluster", cluster)]] <- cov_matrix
  }
  S_pooled <- Reduce("+", lapply(seq_along(cov_matrices), function(i) {
    (group_sizes[i] - 1) * cov_matrices[[i]]
  })) / (sum(group_sizes) - length(unique(res$cluster)))
  
  S_pooled_inv <- solve(S_pooled)
  # Calcul des paramètres w et w0
  gamma <- matrix(0, ncol=G, nrow = ncol(X_tild))
  for (g in 2:G) {
    mu_g <- mu_hat[g, ]
    w_g <- S_pooled_inv %*% mu_g
    w0_g <- -0.5 * t(mu_g) %*% S_pooled_inv %*% mu_g + log(1 / G)  # Assuming equal priors (1/G)
    gamma[,g] <- c(w0_g, w_g)
  }
  
  # Moyenne des lois de y_i
  m_hat <- tapply(y, res$cluster, mean)
  # Variance des loi de y_i
  sigma_hat <- tapply(y, res$cluster, var)
  
  return(list(mu_hat = mu_hat, m_hat = m_hat, gamma = gamma, sigma_h=sigma_hat))
}

init_param_stand <- function(X_tild, y, G) {
  res <- kmeans(X_tild[,-1], centers = G, nstart = 25)
  
  # Moyenne des groupes 
  mu_hat <- res$centers
  
  # Spooled
  cov_matrices <- list()
  group_sizes <- res$size
  for (cluster in unique(res$cluster)) {
    cluster_data <- X_tild[res$cluster == cluster, -1]
    cov_matrix <- cov(cluster_data)
    cov_matrices[[paste("Cluster", cluster)]] <- cov_matrix
  }
  S_pooled <- Reduce("+", lapply(seq_along(cov_matrices), function(i) {
    (group_sizes[i] - 1) * cov_matrices[[i]]
  })) / (sum(group_sizes) - length(unique(res$cluster)))
  
  S_pooled_inv <- solve(S_pooled)
  # Calcul des paramètres w et w0
  gamma <- matrix(0, ncol=G, nrow = ncol(X_tild))
  for (g in 2:G) {
    mu_g <- mu_hat[g, ]
    w_g <- S_pooled_inv %*% mu_g
    w0_g <- -0.5 * t(mu_g) %*% S_pooled_inv %*% mu_g + log(1 / G)
    gamma[,g] <- c(w0_g, w_g)
  }
  
  # Beta_i des lois de y_i
  beta_hat <- matrix(0, nrow = ncol(X_tild), ncol = G)
  for (g in 1:G) {
    group_indices <- which(res$cluster == g)
    model <- lm(y[group_indices] ~ X_tild[group_indices, ] - 1)
    beta_hat[, g] <- coef(model)
  }
  
  # Variance des loi de y_i
  sigma_hat <- tapply(y, res$cluster, var)
  
  return(list(mu = mu_hat, beta = beta_hat, gamma = gamma, sigma=sigma_hat))
}
