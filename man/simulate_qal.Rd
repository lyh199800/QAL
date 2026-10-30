\name{simulate_qal}
\alias{simulate_qal}
\title{Generate an Example from the Supplied QAL Simulation Design}
\usage{simulate_qal(n = 500L, p = 10L, heterogeneous = TRUE)}
\arguments{
\item{n}{Integer sample size at least 20.}
\item{p}{Integer number of covariates at least 6.}
\item{heterogeneous}{Logical; include the treatment-dependent error term when TRUE.}
}
\description{Uses independent uniform covariates, a logistic binary treatment model and a covariate-dependent normal error, following the local revised simulation scripts. Default dimension is smaller for illustration. Controls no random seed; call set.seed explicitly. This helper does not assert an analytic marginal QTE.}
\value{A data frame with X1 through Xp, treatment A and outcome Y.}
\examples{set.seed(2026)
head(simulate_qal())}
