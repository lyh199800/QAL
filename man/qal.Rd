\name{qal}
\alias{qal}
\alias{print.qal}
\alias{coef.qal}
\alias{predict.qal}
\title{Quantile Adaptive Lasso for Binary Treatment}
\usage{
qal(data, covariates, treatment, outcome, tau = .25, known = NULL,
  lambda_exponents = c(-10, -5, -2, -1, -.75, -.5, -.25, .25, .49),
  gamma_convergence = 2, beta_floor = 1e-8,
  selection_digits = 5, ps_clip = NULL)
\method{print}{qal}(x, ...)
\method{coef}{qal}(object, ...)
\method{predict}{qal}(object, newdata = c(0, 1), ...)
}
\arguments{
\item{data}{Data frame containing complete finite numeric analysis columns.}
\item{covariates}{Distinct candidate covariate column names.}
\item{treatment, outcome}{Single treatment and outcome column names. Treatment is coded 0/1.}
\item{tau}{One quantile strictly between zero and one. Fit separately for each target quantile.}
\item{known}{Optional covariate names forced into all adjustment sets, disjoint from candidates.}
\item{lambda_exponents}{Distinct finite exponents e defining lambda = n raised to e.}
\item{gamma_convergence}{Convergence constant c defining gamma = 2(c-e+1).}
\item{beta_floor}{Positive floor on absolute full outcome coefficients for numerical stability.}
\item{selection_digits}{Number of decimal places for rounding treatment coefficients before selection.}
\item{ps_clip}{Optional binary propensity clipping threshold in (0, 0.5); NULL disables clipping.}
\item{x, object}{A fitted qal object.}
\item{newdata}{Numeric vector of treatment values 0 or 1.}
\item{...}{Reserved additional arguments.}
}
\description{Constructs adaptive penalty weights at the requested outcome quantile, tunes propensity covariate selection by balance, and fits a weighted marginal quantile model.}
\details{
Standardizes all covariates using sample standard deviations. An unpenalized full quantile regression of outcome on treatment and covariates provides adaptive penalty weights. Binary treatment selection uses the LQA adaptive lasso. Selected variables enter an unpenalized logistic refit. Inverse probability weights are 1/PS for treated and 1/(1-PS) for controls. An empty set uses unit weights, matching the supplied research implementation.
Tuning minimizes the sum of absolute weighted mean differences in standardized covariates multiplied by absolute full quantile outcome coefficients. The first minimum wins. Failed candidates are recorded and excluded with a warning; all failures cause an error.
The final saturated binary treatment quantile model estimates a difference in marginal potential-outcome quantiles under causal identification assumptions. It does not estimate the quantile of individual treatment effects. Causal interpretation requires consistency, unconfoundedness and positivity.
The full outcome model requires a full-rank design with n greater than covariate count plus two. Numeric indicator coding and any necessary screening must be performed explicitly before fitting.
No standard error or confidence interval is supplied: generic quantile regression summaries do not propagate selection and estimated-weight uncertainty. Inference would require a separately justified procedure, such as resampling the complete analysis pipeline and accounting for failed fits.
}
\value{A qal object with tau, estimate (QTE), quantiles (control and treated), selected variable names, lambda and lambda_exponent, weights, propensity when available, balance diagnostics, tuning scores and errors, candidate path, weighted quantile fit, full outcome_model, centers, scales, ps_clip and effective_sample_size. Coef returns the marginal quantile intercept and treatment coefficient; predict returns the fitted group quantiles.}
\references{Liu Y, Wei K, Huang C, Yu Y, Qin G (2025). Quantile adaptive lasso: variable selection for quantile treatment effect estimation. Journal of Statistical Computation and Simulation 95(2), 239--257. \doi{10.1080/00949655.2024.2424348}. Published online 2024.}
\examples{
set.seed(2026)
d <- simulate_qal(n = 300, p = 6)
f <- qal(d, paste0("X", 1:6), "A", "Y", tau = .25)
f
f$balance
predict(f)
}
