# ----
euclidean_distance <- function(v1, v2) {
  sqrt(sum((v1 - v2)^2))
}

# ----
one_hot_errormachine <- function(Z,size=NULL){
  n <- length(Z)
  if(is.null(size)) K <- length(unique(Z)) else K = size
  mat <- matrix(.Machine$double.eps,n,K)
  mat[cbind(1:n,Z)] <- 1-.Machine$double.eps
  return(mat)
}

# ----
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

# init ----
init_param <- function(X, Q, L) {
  #Tau_mat <- matrix(1/Q, dim(X)[1], Q)
  #Nu_mat <- matrix(1/L, dim(X)[2], L)
  
  #Tau_mat <- matrix(runif(dim(X)[1]*Q), dim(X)[1], Q)
  #Tau_mat <- Tau_mat / apply(Tau_mat,1, sum)
  
  tmp1 <- kmeans(X, Q, nstart = 100)
  cluster1 <- tmp1$cluster
  Tau_mat <- one_hot_errormachine(cluster1, Q)
  alpha <- unname(table(tmp1$cluster)/dim(X)[1])
  
  tmp2 <- kmeans(t(X), L, , nstart = 100)
  beta <- unname(table(tmp2$cluster)/dim(X)[2])
  cluster2 <- tmp2$cluster
  Nu_mat <- one_hot_errormachine(cluster2, L)
  
  
  theta <- array(NA, dim = c(Q,L,2))
  for(i in 1:Q) {
    for (j in 1:L) {
      group <- X[which(cluster1 == i), which(cluster2 == j)]
      mu <- mean(group)
      sig <- sd(group)
      
      theta[i,j,] <- c(mu, sig)
    }
  }
  
  return(list(Tau_mat=Tau_mat, Nu_mat=Nu_mat, alpha=alpha, beta=beta, theta=theta))
}

# init ----
co_init_param <- function(X, Q, L) {
  tmp1 <- kmeans(X, Q, nstart = 500)
  cluster1 <- tmp1$cluster
  Tau_mat <- one_hot_errormachine(cluster1, Q)
  alpha <- unname(table(tmp1$cluster)/dim(X)[1])
  
  clusters_2 <- matrix(NA, Q, dim(X)[2])
  betas <- matrix(NA, Q, L)
  Nu_tensor <- array(0, dim=c(ncol(X), L, Q))
  for(k in 1:Q) {
    tmp2 <- kmeans(t(X[cluster1==k,]), L, nstart=500)
    cluster2 <- tmp2$cluster
    clusters_2[k,] <- cluster2
    betas[k,] <- unname(table(cluster2)/dim(X)[2])
    Nu_tensor[,,k] <- one_hot_errormachine(cluster2, L) 
  }
  
  theta <- array(NA, dim = c(Q,L,2))
  for(i in 1:Q) {
    for (j in 1:L) {
      group <- X[which(cluster1 == i), which(clusters_2[i,] == j)]
      mu <- mean(group)
      sig <- sd(group)
      theta[i,j,] <- c(mu, sig)
    }
  }
  
  return(list(Tau_mat=Tau_mat, Nu_tensor=Nu_tensor, alpha=alpha, betas=betas, theta=theta))
}