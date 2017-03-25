# Compute robust (and clustered) standard errors for speedglm

speedglm.cluster <- function(fit, y, X, cluster = NULL){
  require(Matrix)
  
  if(class(fit) != "speedglm"){
    stop("You must provide a speedglm object.")
  }
  
  coefs <- fit$coefficients
  X <- X[, !is.na(coefs)]
  coefs <- coefs[!is.na(coefs)]
  n <- fit$n
  k <- length(coefs)
  eta <- as.matrix(X %*% coefs)
  mu <- fit$family$linkinv(eta)
  res <- y - mu
  
  W1 <- Diagonal(mu, n = n)
  bread <- solve(t(X) %*% W1 %*% X)
  
  if(is.null(cluster)){
    W2 <- Diagonal(res^2, n = n)
    meat <- t(X) %*% W2 %*% X
    
  } else {
    m <- length(unique(cluster))
    dfc <- (m/(m-1))*((n-1)/(n-k))
    
    meat <- matrix(0, nrow = k, ncol = k)
    
    for(i in unique(cluster)){
      X_tmp <- as.matrix(X[cluster == i, ])
      res_tmp <- res[cluster == i]
      W2_tmp <- res_tmp %*% t(res_tmp)
      meat <- meat + (t(X_tmp) %*% W2_tmp %*% X_tmp)
    }
    
  }
  
  vcov <- bread %*% meat %*% bread
  se <- sqrt(diag(vcov))
  z_value <- coefs/se

  results <- data.frame(Estimate = coefs,
                        Std_Error = se,
                        z_value = coefs/se,
                        p_value = (1 - pnorm(abs(z_value))))

  results
}
