# Quantile adaptive lasso: Liu et al., doi:10.1080/00949655.2024.2424348.
qal <- function(data, covariates, treatment, outcome, tau = .25,
                known = NULL,
                lambda_exponents = c(-10, -5, -2, -1, -.75, -.5, -.25, .25, .49),
                gamma_convergence = 2, beta_floor = 1e-8,
                selection_digits = 5, ps_clip = NULL) {
  call <- match.call()
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!is.character(covariates) || !length(covariates) || anyDuplicated(covariates))
    stop("covariates must contain distinct column names.")
  if (!is.null(known) && (!is.character(known) || anyDuplicated(known)))
    stop("known must contain distinct column names.")
  if (!is.character(treatment) || length(treatment) != 1L ||
      !is.character(outcome) || length(outcome) != 1L)
    stop("treatment and outcome must each name one column.")
  cols <- c(known, covariates, treatment, outcome)
  if (anyDuplicated(cols) || anyDuplicated(names(data)) || !all(cols %in% names(data)))
    stop("Analysis column names must exist, be unique and non-overlapping.")
  if (!all(vapply(data[cols], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[cols]))))
    stop("Analysis columns must be numeric and finite; handle missing values explicitly.")
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(tau) || tau <= 0 || tau >= 1) stop("tau must be one number strictly between 0 and 1.")
  if (!is.numeric(lambda_exponents) || !length(lambda_exponents) ||
      any(!is.finite(lambda_exponents)) || anyDuplicated(lambda_exponents))
    stop("lambda_exponents must contain distinct finite numbers.")
  if (!scalar(gamma_convergence) || any(2 * (gamma_convergence - lambda_exponents + 1) <= 0))
    stop("gamma_convergence must yield positive adaptive exponents.")
  if (!scalar(beta_floor) || beta_floor <= 0) stop("beta_floor must be positive.")
  if (!scalar(selection_digits) || selection_digits < 0 || selection_digits > 12 ||
      selection_digits %% 1 != 0) stop("selection_digits must be an integer from 0 to 12.")
  if (!is.null(ps_clip) && (!scalar(ps_clip) || ps_clip <= 0 || ps_clip >= .5))
    stop("ps_clip must be NULL or a number strictly between 0 and 0.5.")
  xn <- c(known, covariates); X <- as.matrix(data[xn])
  A <- data[[treatment]]; Y <- data[[outcome]]; n <- nrow(X)
  if (!setequal(unique(A), c(0, 1))) stop("Treatment must contain both 0 and 1.")
  if (n <= ncol(X) + 2L) stop("Full quantile model requires n > number of covariates + 2; screen first.")
  scales <- apply(X, 2, stats::sd); centers <- colMeans(X)
  if (any(scales <= 0)) stop("Remove constant covariates.")
  X <- scale(X, center = centers, scale = scales)
  if (qr(cbind(1, A, X))$rank < ncol(X) + 2L) stop("Remove collinear covariates or treatment predictors.")
  internal <- paste0("x", seq_len(ncol(X)))
  d <- as.data.frame(X); names(d) <- internal; d$a <- A; d$y <- Y
  outcome_model <- quantreg::rq(stats::reformulate(c("a", internal), "y"), data = d, tau = tau)
  beta <- stats::coef(outcome_model)[internal]
  if (any(!is.finite(beta))) stop("Non-finite full quantile coefficients.")
  names(beta) <- xn
  M <- stats::model.matrix(stats::reformulate(internal), d)
  paths <- vector("list", length(lambda_exponents)); rows <- paths
  for (i in seq_along(paths)) {
    exponent <- lambda_exponents[i]; gamma <- 2 * (gamma_convergence - exponent + 1)
    lambda <- n^exponent
    ans <- tryCatch({
      pw <- exp(-gamma * log(pmax(abs(beta), beta_floor)))
      if (length(known)) pw[seq_along(known)] <- 0
      if (any(!is.finite(pw)) || !is.finite(lambda) || lambda <= 0) stop("Invalid adaptive penalty.")
      penalty <- lqa::adaptive.lasso(lambda = lambda, al.weights = pw)
      sel <- lqa::lqa.default(x = M, y = A, penalty = penalty, family = stats::binomial())
      b <- stats::coef(sel)[internal]
      if (any(!is.finite(b))) stop("Non-finite treatment coefficients.")
      keep <- sort(unique(c(seq_along(known), which(round(b, selection_digits) != 0))))
      if (length(keep)) {
        psfit <- stats::glm(stats::reformulate(internal[keep], "a"), data = d, family = stats::binomial())
        if (!isTRUE(psfit$converged)) stop("Propensity refit did not converge.")
        ps <- stats::predict(psfit, type = "response")
        if (!is.null(ps_clip)) ps <- pmin(1 - ps_clip, pmax(ps_clip, ps))
        if (any(!is.finite(ps)) || any(ps <= 0 | ps >= 1)) stop("Boundary propensity scores; inspect overlap.")
        w <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
      } else { w <- rep(1, n); ps <- NULL }
      if (any(!is.finite(w)) || any(w <= 0)) stop("Invalid inverse probability weights.")
      balance <- .qal_balance(X, A, w)
      score <- sum(abs(beta) * balance)
      if (!is.finite(score)) stop("Non-finite balance score.")
      list(selected = xn[keep], weights = w, propensity = ps, balance = balance,
           score = score, treatment_coef = stats::setNames(b, xn))
    }, error = function(e) list(error = conditionMessage(e)))
    paths[[i]] <- ans
    rows[[i]] <- data.frame(exponent = exponent, lambda = lambda, gamma = gamma,
      score = if (is.null(ans$error)) ans$score else Inf,
      n_selected = if (is.null(ans$error)) length(ans$selected) else NA_integer_,
      error = if (is.null(ans$error)) "" else ans$error)
  }
  tuning <- do.call(rbind, rows); good <- is.finite(tuning$score)
  if (!any(good)) stop("All tuning candidates failed: ", paste(unique(tuning$error), collapse = "; "))
  if (any(!good)) warning(sum(!good), " tuning candidates failed; inspect $tuning.")
  best <- which.min(tuning$score); chosen <- paths[[best]]
  fit <- quantreg::rq(y ~ a, data = d, tau = tau, weights = chosen$weights)
  cf <- stats::coef(fit)
  if (any(!is.finite(cf))) stop("Non-finite final weighted quantile coefficients.")
  result <- list(call = call, tau = tau, estimate = unname(cf["a"]),
    quantiles = c(control = unname(cf[1]), treated = unname(sum(cf))),
    selected = chosen$selected, lambda_exponent = tuning$exponent[best],
    lambda = tuning$lambda[best], weights = chosen$weights, propensity = chosen$propensity,
    balance = data.frame(covariate = xn,
      unweighted = as.numeric(.qal_balance(X, A, rep(1, n))),
      weighted = as.numeric(chosen$balance), outcome_coefficient = as.numeric(beta)),
    tuning = tuning, path = paths, fit = fit, outcome_model = outcome_model,
    centers = centers, scales = scales, ps_clip = ps_clip,
    effective_sample_size = sum(chosen$weights)^2 / sum(chosen$weights^2))
  class(result) <- "qal"; result
}

.qal_balance <- function(X, A, w) {
  apply(X, 2, function(x) abs(stats::weighted.mean(x[A == 1], w[A == 1]) -
    stats::weighted.mean(x[A == 0], w[A == 0])))
}
print.qal <- function(x, ...) {
  cat("QAL: tau =", x$tau, "\nQuantile treatment effect:", format(x$estimate, digits = 5), "\n")
  cat("Selected:", if (length(x$selected)) paste(x$selected, collapse = ", ") else "(none)", "\n")
  cat("Lambda exponent:", x$lambda_exponent, "  Effective sample size:", round(x$effective_sample_size, 1), "\n")
  invisible(x)
}
coef.qal <- function(object, ...) stats::coef(object$fit)
predict.qal <- function(object, newdata = c(0, 1), ...) {
  if (!is.numeric(newdata) || any(!is.finite(newdata)) || any(!newdata %in% c(0, 1)))
    stop("newdata must contain binary treatment values 0 or 1.")
  b <- stats::coef(object$fit); unname(b[1] + b[2] * newdata)
}
