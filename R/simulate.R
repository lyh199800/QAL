simulate_qal <- function(n = 500L, p = 10L, heterogeneous = TRUE) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 20 || n %% 1 != 0 ||
      !is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 6 || p %% 1 != 0)
    stop("Use integer n >= 20 and p >= 6.")
  if (!is.logical(heterogeneous) || length(heterogeneous) != 1L || is.na(heterogeneous))
    stop("heterogeneous must be TRUE or FALSE.")
  X <- matrix(stats::runif(n * p), n, p)
  v <- stats::rnorm(n)
  c0 <- -1 / stats::qnorm(.25)
  A <- stats::rbinom(n, 1, stats::plogis(rowSums(X[, 1:5, drop = FALSE])))
  z <- if (heterogeneous) .1 else 0
  Y <- 1 + 2 * A + rowSums(X[, 1:6, drop = FALSE]) +
    (c0 * (X[, 2] + X[, 3] + X[, 5]) + z * A) * v
  d <- as.data.frame(X); names(d) <- paste0("X", seq_len(p)); d$A <- A; d$Y <- Y
  d
}
