gravity_ppml <- function(y, x, fixed_effects, data,
                         offset = NULL,
                         subset = NULL,
                         robust = TRUE,
                         cluster = NULL,
                         tol = 1e-10) {
  if (!is.null(subset)) data <- data[subset, ]

  offset2 <- offset
  if (is.null(offset)) offset <- rep(0, nrow(data))

  trade <- data[[y]]

  if (min(trade) < 0) {
    stop("y should be greater or equals to zero.")
  }

  data[, fixed_effects] <- lapply(data[, fixed_effects], as.factor)

  max_trade <- max(trade)
  trade <- trade / max_trade

  mu <- (trade + 0.5) / 2
  eta <- log(mu) - offset
  z <- eta + (trade - mu) / mu

  # Formula

  if (!is.null(x)) {
    xvars <- paste0(x, collapse = " + ")
  }

  f <- paste0(
    "z ~ ", ifelse(is.null(x), " -1 ", xvars), " | ",
    paste0(fixed_effects, collapse = " + ")
  )
  f <- as.formula(f)

  dif <- 1
  rss1 <- 1
  while (abs(dif) > tol) {
    reg <- lfe::felm(f,
      data = data,
      weights = mu
    )

    eta <- z - reg$residuals + offset
    mu <- exp(eta)
    z <- (eta - offset) + (trade - mu) / mu

    res <- trade - mu
    rss2 <- sum(res^2)
    dif <- rss2 - rss1
    rss1 <- rss2
    dev <- 2 * max_trade * sum(trade[trade > 0] * log(trade[trade > 0] / mu[trade > 0]))
    # cat(dev, "\n")
    # cat('Sum of Squared Residuals = ', rss2, "\n")
  }

  z <- z + log(max_trade)
  reg <- lfe::felm(f,
    data = data,
    weights = mu
  )

  if (!is.null(x)) {
    z <- data.frame(id = 1:nrow(data))
    for (i in x) {
      fixed_effects_tmp <- paste0(fixed_effects, collapse = " + ")
      f <- as.formula(paste0(
        i, " ~ -1 ",
        ifelse(!is.null(offset2), " + offset ", ""),
        "| ", fixed_effects_tmp, " | 0 | 0"
      ))
      fit.tmp <- lfe::felm(f, data = data, weights = mu)
      z[[i]] <- fit.tmp$residuals
    }
    z <- z[, -1]
    z <- as.matrix(z)

    n <- reg$N
    k <- length(x)
    W1 <- Matrix::Diagonal(mu, n = n)
    bread <- solve(t(z) %*% W1 %*% z)

    cluster_name <- cluster

    res <- trade - mu
    if (robust) {
      if (is.null(cluster)) {
        W2 <- Diagonal((res^2), n = n)
        meat <- t(z) %*% W2 %*% z
      } else {
        cluster <- data[[cluster]]
        m <- length(unique(cluster))
        dfc <- (m / (m - 1)) * ((n - 1) / (n - k))

        meat <- matrix(0, nrow = length(x), ncol = length(x))

        for (i in unique(cluster)) {
          z_tmp <- as.matrix(z[cluster == i, , drop = FALSE])
          res_tmp <- res[cluster == i]
          W2_tmp <- res_tmp %*% t(res_tmp)
          meat <- meat + (t(z_tmp) %*% W2_tmp %*% z_tmp)
        }
      }
    }

    vcov <- if (robust) {
      bread %*% meat %*% bread
    } else {
      bread
    }

    reg$vcv <- vcov
    reg$se <- sqrt(diag(reg$vcv))
    reg$tval <- reg$coefficients / reg$se
    reg$pval <- 1 - pnorm(abs(reg$tval))

    if (robust) {
      reg$vcv <- vcov
      reg$se <- sqrt(diag(reg$vcv))
      reg$tval <- reg$coefficients / reg$se
      reg$pval <- 1 - pnorm(abs(reg$tval))

      reg$robustvcv <- vcov
      reg$rse <- sqrt(diag(reg$robustvcv))
      reg$rtval <- reg$coefficients / reg$rse
      reg$rpval <- 1 - pnorm(abs(reg$rtval))
    }
  }

  x_fixed_effects <- data[, fixed_effects]
  x_fixed_effects$order <- 1:nrow(x_fixed_effects)
  len_fe <- length(fixed_effects)

  for (i in 1:len_fe) {
    fe_tmp <- lfe::getfe(reg)
    fe_tmp <- fe_tmp[fe_tmp$fe == fixed_effects[i], c("idx", "effect")]

    colnames(fe_tmp) <- c(
      fixed_effects[i],
      paste0("fe_", fixed_effects[i])
    )
    x_fixed_effects <- merge(x_fixed_effects,
      fe_tmp,
      by = fixed_effects[i],
      all.x = TRUE
    )
  }
  x_fixed_effects <- x_fixed_effects[
    order(x_fixed_effects$order),
    -(len_fe + 1)
  ]
  x_fixed_effects[, 1:len_fe] <- sapply(x_fixed_effects[, 1:len_fe], as.character)
  reg$fixed.effects <- x_fixed_effects

  x_fixed_effects <- x_fixed_effects[, !names(x_fixed_effects) %in% fixed_effects]
  x_fixed_effects <- apply(x_fixed_effects, 1, sum)

  if (!is.null(x)) {
    x_var <- as.matrix(data[, x])
    beta <- as.matrix(reg$coefficients)
    x_beta <- exp(x_var %*% reg$coefficients + offset + x_fixed_effects)
  } else {
    x_beta <- exp(offset + x_fixed_effects)
  }

  reg$fitted.values <- x_beta

  reg$R2 <- cor(data[[y]], x_beta)^2

  class(reg) <- "gravity.ppml"
  return(reg)
}

summary.gravity.ppml <- function(object) {
  class(object) <- "summary.gravity.ppml"
  object
}

print.summary.gravity.ppml <- function(object) {
  cat("Coefficients: \n")
  results <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `t value` = object$tval,
    `Pr(>|t|)` = object$pval
  )
  results <- as.matrix(results)
  colnames(results) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(results, digits = 4)
}

predict.gravity.ppml <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(offset)) offset <- rep(0, nrow(newdata))

  fixed_effects <- names(object$fe)
  x_fixed_effects <- newdata[, fixed_effects]
  x_fixed_effects$order <- 1:nrow(x_fixed_effects)
  len_fe <- length(fixed_effects)

  for (i in 1:len_fe) {
    fe_tmp <- getfe(object)
    fe_tmp <- fe_tmp[fe_tmp$fe == fixed_effects[i], c("idx", "effect")]

    colnames(fe_tmp) <- c(
      fixed_effects[i],
      paste0("fe_", fixed_effects[i])
    )

    x_fixed_effects <- merge(x_fixed_effects,
      fe_tmp,
      by = fixed_effects[i],
      all.x = TRUE
    )
  }

  x_fixed_effects <- x_fixed_effects[
    order(x_fixed_effects$order),
    -(len_fe + 1)
  ]

  x_fixed_effects[, 1:len_fe] <- sapply(x_fixed_effects[, 1:len_fe], as.character)
  object$fixed.effects <- x_fixed_effects

  x_fixed_effects <- x_fixed_effects[, !names(x_fixed_effects) %in% fixed_effects]
  x_fixed_effects <- apply(x_fixed_effects, 1, sum)

  x <- rownames(object$beta)

  if (!is.null(x)) {
    x_var <- as.matrix(newdata[, x])
    beta <- as.matrix(object$coefficients)
    x_beta <- exp(x_var %*% object$coefficients + offset + x_fixed_effects)
  } else {
    x_beta <- exp(offset + x_fixed_effects)
  }

  x_beta
}