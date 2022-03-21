
compute_gravity_vars <- function(p, P, PI, Y, E, sigma, tc, N) {
  N <- length(p)

  X <- matrix(nrow = N, ncol = N)

  for (j in 1:N) {
    X[, j] <- (Y * E[j] / sum(Y)) * (tc[, j] / (P[j] * PI))^(1 - sigma)
  }

  list(
    p = p,
    Y = Y,
    E = E,
    P = P,
    PI = PI,
    X = X,
    RGDP = Y / P
  )
}

solve_gravity <- function(start, phi, Q, alpha, tc, N) {
  p <- start[1:N]
  P <- start[(N + 1):(2 * N)]
  PI <- start[(2 * N + 1):(3 * N)]

  Y <- p * Q
  E <- phi * p * Q

  PI_new <- PI
  P_new <- P

  for (i in 1:N) {
    PI_new[i] <- sum((tc[i, ] / P)^(1 - sigma) * E / sum(Y))^(1 / (1 - sigma))
  }

  for (j in 1:N) {    
    P_new[j] <- sum((tc[, j] / PI)^(1 - sigma) * Y / sum(Y))^(1 / (1 - sigma))
  }

  X <- compute_gravity_vars(p, P, PI, Y, E, sigma, tc, N)$X

  res <- c(
    p - (Y / sum(Y))^(1 / (1 - sigma)) * 1 / (alpha * PI),
    P[1] - 1,
    P[-1] - P_new[-1],
    PI - PI_new
  )

  res
}

solve_mrts <- function(start, N, tc, sigma, E, Y) {
  P <- start[1:N]
  PI <- start[(N + 1):(2 * N)]

  PI_new <- PI
  P_new <- P

  for (i in 1:N) {
    PI_new[i] <- sum((tc[i, ] / P)^(1 - sigma) * E / sum(Y))^(1 / (1 - sigma))
  }

  for (j in 1:N) {
    P_new[j] <- sum((tc[, j] / PI)^(1 - sigma) * Y / sum(Y))^(1 / (1 - sigma))
  }

  c(
    P[1] - 1,
    P[2:N] - P_new[-1],
    PI - PI_new
  )
}

ge_gravity <- function(phi, Q, alpha, tc, N, full_endowment = TRUE) {
  PI <- (Y / sum(Y))^(1 / (1 - sigma)) * 1 / alpha
  P <- PI
  for (j in 1:N) {
    P[j] <- sum((tc[, j] / PI)^(1 - sigma) * Y / sum(Y))^(1 / (1 - sigma))
  }


  if (full_endowment) {
    solution <- nleqslv::nleqslv(
      x = c(rep(1, N), P, PI),
      fn = solve_gravity,
      phi = phi,
      Q = Q,
      alpha = alpha,
      tc = tc,
      N = N
    )
    if (solution$message != "Function criterion near zero") {
      stop("Convergence problem!")
    }
    p <- solution$x[1:N]
    P <- solution$x[(N + 1):(2 * N)]
    PI <- solution$x[(2 * N + 1):(3 * N)]
    Y <- p * Q
    E <- phi * Y
  } else {
    p <- rep(1, N)
    Y <- p * Q
    E <- phi * Y

    sol_mrts <- nleqslv::nleqslv(
      x = c(P, PI),
      fn = solve_mrts,
      N = N,
      tc = tc,
      sigma = sigma,
      E = E,
      Y = Y
    )
    if (sol_mrts$message != "Function criterion near zero") {
      stop("Convergence problem!")
    }
    P <- sol_mrts$x[1:N]
    PI <- sol_mrts$x[(N + 1):(2 * N)]
  }

  gravity_vars <- compute_gravity_vars(p, P, PI, Y, E, sigma, tc, N)

  gravity_vars
}