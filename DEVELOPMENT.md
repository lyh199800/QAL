# Source mapping and validation

Sources inspected: Rcode_Revise/QAL_function.R, QAL_function_gamma5.R,
Settings_Q_gamma2_qtediff_25th_500.R and the revised manuscript supplied locally.
The authoritative method reference is doi:10.1080/00949655.2024.2424348.

## Preserved algorithm

Sample-standardized covariates; full outcome quantile regression at tau;
adaptive LQA binomial treatment model; lambda=n^e; gamma=2(c-e+1);
five-decimal coefficient selection; logistic propensity refit; binary IPW;
wAMD with absolute quantile outcome coefficients; first minimum on the grid;
weighted quantile regression of outcome on binary treatment. Empty sets use
unit weights. Default c=2. The original gamma5 file also omits -10 from its
exponent grid; match both c=5 and that shorter grid to reproduce that version.

## Explicit changes

- Explicit data and tau replace global Y, tau and var.list dependencies.
- Misleading GOAL/ATE function and result names are replaced with qal/QTE names.
- Direct LQA matrix interface avoids the archived formula contrasts warning.
- beta_floor=1e-8 avoids infinite penalties for zero outcome coefficients.
- Input checks reject non-finite data, constant columns, invalid binary coding,
  insufficient observations and rank deficiency. Candidates that fail are
  logged and excluded with a warning; if all fail the fit stops.
- Optional known covariates extend the old unused covar argument: always
  included with zero penalty. They also enter balance tuning.
- Propensity clipping is explicit and off by default; enabling it changes weights.
- Survey normalization in ATE_est was a scalar divided by itself (equal to 1).
  It is removed. The final rq receives the same IPW; survey is not a dependency.
- Standard errors and intervals are not invented. Users must justify a full
  selection-and-weight-aware inferential procedure separately.
- Simulator follows the revised heteroskedastic example, using p=10 by default;
  it does not claim that a conditional coefficient formula equals a marginal QTE.

## Validation

Delivery validation runs the actual supplied QAL functions in an isolated
environment, removing library attachments and explicitly supplying required
workspace variables. At tau=.25 and .75, and c=2 and 5, it compares selected
covariates, selected tuning exponent, IPW and final QTE against the new package.
The original input files are not modified or redistributed.
Package tests additionally check balance, final quantile regression, shift
invariance, column renaming, known-covariate inclusion and invalid inputs.
R CMD check uses --no-manual; Windows/Linux checks await GitHub Actions.
