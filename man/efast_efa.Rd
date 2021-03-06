% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efast.R
\name{efast_efa}
\alias{efast_efa}
\title{Estimate an EFA model in lavaan}
\usage{
efast_efa(
  data,
  M,
  sample.nobs = NULL,
  auto.fix.first = FALSE,
  auto.var = TRUE,
  auto.efa = TRUE,
  information = "observed",
  std.ov = TRUE,
  ...
)
}
\arguments{
\item{data}{<data.frame> the dataset or <matrix> covariance matrix}

\item{M}{<numeric> How many factors, minimum 2}

\item{sample.nobs}{<numeric> sample size (if data = covmat, see lavaan)}

\item{auto.fix.first}{<bool> see lavaan}

\item{auto.var}{<bool> see lavaan}

\item{auto.efa}{<bool> see lavaan}

\item{information}{<character> see lavaan}

\item{std.ov}{<bool> see lavaan}

\item{...}{other arguments passed to lavaan}
}
\description{
This function estimates efa models in the same way that efast models are
estimated: using lavaan.
}
\details{
The constrained model constrains the residual covariance to be
  equal across the different ROIs.
}
\examples{
\dontrun{
# create a test dataset
test_data <- simulate_efast()
fit_efa <- efast_efa(simdat, M = 4)
summary(fit_efa)
}

}
