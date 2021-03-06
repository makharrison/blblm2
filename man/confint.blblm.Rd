% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{confint.blblm}
\alias{confint.blblm}
\title{Confidence interval for blblm coefficients}
\usage{
\method{confint}{blblm}(object, parm = NULL, level = 0.95, ...)
}
\arguments{
\item{object}{A fitted blblm model}

\item{parm}{The number of parameters in the fitted blblm model. Default is set to NULL.}

\item{level}{A numeric value indicating the desired level of confidence interval.}

\item{...}{Additional parameters}
}
\value{
A confidence interval for the blblm coefficients.
}
\description{
Compute confidence interval for the coefficients of a blblm model.
}
