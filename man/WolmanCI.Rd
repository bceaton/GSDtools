% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSDtools_functions.r
\name{WolmanCI}
\alias{WolmanCI}
\title{Calculate the percentile sizes and grain size confidence intervals for a sample}
\usage{
WolmanCI(cfd, n, P = seq(5, 95, 5), equaltail = T, alpha = 0.05)
}
\arguments{
\item{cfd}{is a data frame providing a list grain sizes in the first variable
and the corresponding cumulative proportion finer in the second. The
grain sizes should be recorded in mm, and the proportion finer in [0,1].}

\item{n}{is the total number of observations upon which the cumulative
frequency distribution in \code{cfd} is based}

\item{P}{numeric vector of percentiles to estimate. The default is to estimate
the 5th to the 95th percentile, in increments of 5.}

\item{equaltail}{is a logical variable that determines whether the calculations
of the confidence interval are based on an approximation with equal areas in
each tail (the default), or based on the exact binomial solution with a coverage
of at least 95\%}

\item{alpha}{the desired confidence level for which to calculate a
confidence interval in [0,1].}
}
\description{
\code{WolmanCI} is a function that uses cumulative frequency distribution
data for Wolman or grid-by-number samples of bed surface texture to
estimate the value of user-specified percentiles. The function also uses
the binomial distribution to estimate the confidence intervals corresponding
to a user-specified confidence level, based on the number of grain size
measurements that were used to construct the cumulative frequency
distribution.
}
\details{
The function returns a data frame listing the estimate of each percentile,
the upper limit, and the lower limit.
}
