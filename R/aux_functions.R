#
gravity_ppml <- function(y, x, data, fixed_effects = NULL, cluster = NULL,
                         reference = NULL, subset = NULL, robust = TRUE){
  require(Matrix)
  require(speedglm)
  
  if(!is.null(subset)) data <- data[subset,]
  
  xvars <- c(x, fixed_effects)
  xvars <- paste0(xvars, collapse = " + ")
  x_formula <- as.formula(paste0("~ -1 + ", xvars))
  X <- sparse.model.matrix(x_formula, data)
  y_data <- data[[y]]
  
  if(!is.null(reference)){
    r <- colnames(X)
    idx <- which(!grepl(reference, r, fixed = TRUE))
    X <- X[, idx]
  }
  
  fit <- speedglm.wfit(y      = y_data,
                       X      = X,
                       sparse = TRUE,
                       family = quasipoisson(link = log),
                       fitted = TRUE)
  
  # Standard Errors
  
  coefs <- fit$coefficients
  X <- X[, !is.na(coefs)]
  coefs <- coefs[!is.na(coefs)]
  n <- fit$n
  k <- length(coefs)
  eta <- as.matrix(X %*% coefs)
  mu <- fit$family$linkinv(eta)
  res <- y_data - mu
  
  W1 <- Diagonal(mu, n = n)
  bread <- solve(t(X) %*% W1 %*% X)
  
  if(robust){
    if(is.null(cluster)){
      W2 <- Diagonal(res^2, n = n)
      meat <- t(X) %*% W2 %*% X
      
    } else {
      cluster <- data[[cluster]]
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
    
    fit$vcov <- bread %*% meat %*% bread
    fit$robust <- robust
    fit$x <- x
    
  }
  
  return(fit)
}

summary.ppml <- function(fit){
  if(fit$robust){
    vcov <- fit$vcov
    x <- fit$x
    sd <- rep(NA_real_, length(fit$coefficients))
    sd[!is.na(fit$coefficients)] <- sqrt(diag(vcov))
    
    s <- summary(fit)
    idx_x <- which(names(fit$coefficients) %in% x)
    s$coefficients <- s$coefficients[idx_x,]
    s$coefficients$`Std. Error` <- sd[idx_x]
    s$coefficients$`t value` <- s$coefficients$Estimate/s$coefficients$`Std. Error`
    s$coefficients$`Pr(>|t|)` <- 1 - pnorm(abs(s$coefficients$`t value`))
    
    s
  } else{
    summary(fit)
  }
}

#
gravity_ppml2 <- function(y, x, data, fixed_effects = NULL, cluster = NULL,
                         reference = NULL, subset = NULL, robust = TRUE, ...){
  require(Matrix)
  require(speedglm)
  require(lfe)
  
  if(!is.null(subset)) data <- data[subset,]
  
  xvars <- c(x, fixed_effects)
  xvars <- paste0(xvars, collapse = " + ")
  x_formula <- as.formula(paste0("~ -1 + ", xvars))
  X <- sparse.model.matrix(x_formula, data)
  y_data <- data[[y]]
  
  if(!is.null(reference)){
    r <- colnames(X)
    idx <- which(!grepl(reference, r, fixed = TRUE))
    X <- X[, idx]
  }
  
  fit <- speedglm.wfit(y      = y_data,
                       X      = X,
                       sparse = TRUE,
                       family = quasipoisson(link = log),
                       fitted = TRUE, ...)
  
  # Standard Errors
  coefs <- fit$coefficients
  X <- X[, !is.na(coefs)]
  coefs <- coefs[!is.na(coefs)]
  n <- fit$n
  k <- length(coefs)
  eta <- as.matrix(X %*% coefs)
  mu <- fit$family$linkinv(eta)
  res <- y_data - mu
  
  z <- data.frame(id = 1:nrow(data))
  for(i in x){
    fixed_effects_tmp <- paste0(fixed_effects, collapse = " + ")
    f <- as.formula(paste0(i, " ~ -1 | ", fixed_effects_tmp, " | 0 | 0"))
    fit.tmp <- felm(f, data = data, weights = mu)
    z[[i]] <- fit.tmp$residuals
  }
  z <- z[,-1]
  
  W1 <- Diagonal(mu, n = n)
  bread <- solve(t(z) %*% W1 %*% z)
  
  if(robust){
    if(is.null(cluster)){
      W2 <- Diagonal((res^2), n = n)
      meat <- t(z) %*% W2 %*% z
      
    } else {
      cluster <- data[[cluster]]
      m <- length(unique(cluster))
      dfc <- (m/(m-1))*((n-1)/(n-k))
      
      meat <- matrix(0, nrow = length(x), ncol = length(x))
      
      for(i in unique(cluster)){
        z_tmp <- as.matrix(z[cluster == i, ])
        res_tmp <- res[cluster == i]
        W2_tmp <- res_tmp %*% t(res_tmp)
        meat <- meat + (t(z_tmp) %*% W2_tmp %*% z_tmp)
      }
    }
    
    fit$vcov <- bread %*% meat %*% bread
    fit$robust <- robust
    fit$x <- x
    
  }
  
  return(fit)
}

