#' Generate a data frame containing cumulative frequency distribution data
#'
#' \code{MakeCFD} is a function that creates a data frame containing the
#' cumulative grain size distribution information for a sample of n observations
#' of b axis diameter (in mm).
#'
#' The function returns a data frame listing upper bounds of each size class
#' containing data, as well as the cumulative proportion of the observations
#' that fall below the given grain size. The output of this function is the
#' required input for other functions in this package. In order to retain
#' information about the sample size, there is the option to output the total
#' counts for each size class, as well.
#'
#' @param obs vector containing b-axis measurements of gravel bed surface
#' @param increment an optional parameter to control the size of the grain size
#' classes used for the analysis (1.0 = 1.0 phi intervals, 0.5 = 1/2 phi
#' intervals, 0.25 = 1/4 phi intervals, etc.). Default value is 0.5.
#' @param count optional flag to record the number of observations falling in
#' each size class as an additional variable in the output data frame
#' @param plot optional flag to produce a graph of the resulting grain size
#' distribution
#' @export

MakeCFD = function(obs, increment = 0.5, count = FALSE, plot = FALSE){
  n = length(obs)
  max.phi = ceiling(log2(max(obs)))
  min.phi = floor(log2(min(obs)))
  size = 2^seq(from = min.phi, to = max.phi, by = increment)
  results = hist(obs,
                 breaks = size,
                 plot = FALSE)
  cfd = data.frame(size)
  cfd$probs = c(0, cumsum(results$counts))/sum(results$counts)

  if(count == TRUE){
    cfd$counts = c(0,results$counts)
  }

  if(plot == TRUE){
    plot(cfd[[1]], cfd[[2]],
         type = "b",
         col = "blue",
         log = "x",
         xlab = "grain size (mm)",
         ylab = "cum. prop. finer")
  }

  return(cfd)
}

#' Calculate the confidence interval for a quantile 'p'
#'
#' \code{QuantBD} uses the  binomial distribution to compute a
#' near-symmetric distribution-free confidence interval for a quantile 'p'.
#' Returns indices for the order statistics, along with coverage probability.
#'
#'The function returns the indices representing a percentile confidence interval
#'that contains the population quantile of interest, assuming a given
#'confidence level. The script also returns the probability coverage
#'for the confidence interval. Finally, the function returns an approximation
#'of the confidence interval that has equal areas in the tails of the distribution
#'
#'@param n the sample size
#'@param p the desired quantile to be estimated from the sample in the range [0,1]
#'@param alpha the confidence level (default value is 0.05)
#'@export

QuantBD <- function(n, p, alpha = 0.05) {
  u <- qbinom(1 - alpha/2, n, p) + (-2:2) + 1
  l <- qbinom(alpha/2, n, p) + (-2:2)
  u[u > n] <- Inf
  l[l < 0] <- -Inf
  p_c <- outer(l, u, function(a, b) pbinom(b-1, n, p) - pbinom(a-1, n, p))
  if (max(p_c) < 1 - alpha) {
    i <- which(p_c == max(p_c))
  } else {
    i <- which(p_c == min(p_c[p_c >= 1 - alpha]))
  }
  i <- i[1]

  u <- rep(u, each = 5)[i]
  l <- rep(l, 5)[i]

  k = 1:n
  pcum = pbinom(k, n , p)

  lu_approx = approx(x = pcum, y = k, xout = c(alpha/2, 1 - alpha/2))$y

  list(interval = c(l, u), coverage = p_c[i], equaltail = lu_approx)
}


#' Calculate the percentile sizes and grain size confidence intervals for a sample
#'
#' \code{WolmanCI} is a function that uses cumulative frequency distribution
#' data for Wolman or grid-by-number samples of bed surface texture to
#' estimate the value of user-specified percentiles. The function also uses
#' the binomial distribution to estimate the confidence intervals corresponding
#' to a user-specified confidence level, based on the number of grain size
#' measurements that were used to construct the cumulative frequency
#' distribution.
#'
#' The function returns a data frame listing the estimate of each percentile,
#' the upper limit, and the lower limit.
#'
#' @param cfd is a data frame providing a list grain sizes in the first variable
#' and the corresponding cumulative proportion finer in the second. The
#' grain sizes should be recorded in mm, and the proportion finer in [0,1].
#' @param n is the total number of observations upon which the cumulative
#' frequency distribution in \code{cfd} is based
#' @param P numeric vector of percentiles to estimate. The default is to estimate
#' the 5th to the 95th percentile, in increments of 5.
#' @param equaltail is a logical variable that determines whether the calculations
#' of the confidence interval are based on an approximation with equal areas in
#' each tail (the default), or based on the exact binomial solution with a coverage
#' of at least 95\%
#' @param alpha  the desired confidence level for which to calculate a
#' confidence interval in [0,1].
#' @export

WolmanCI = function(cfd, n,  P = seq(5, 95, 5), equaltail = T,  alpha = 0.05){
  # use the binomial approach
  probs = P/100  #convert percentiles to probabilities
  p.upper = vector(mode="numeric", length = length(probs))
  p.lower = vector(mode="numeric", length = length(probs))
  if (equaltail == T){
    for(i in seq_along(probs)){
      tmp = QuantBD(n, probs[i], alpha)
      p.upper[i] = tmp$equaltail[2]/n
      p.lower[i] = tmp$equaltail[1]/n
    }
  }else{
    for(i in seq_along(probs)){
      tmp = QuantBD(n, probs[i], alpha)
      p.upper[i] = tmp$interval[2]/n
      p.lower[i] = tmp$interval[1]/n
    }
  }

  # estimate percentiles
  phi = log2(cfd[[1]])
  X = cfd[[2]]
  estimate = 2^approx(x = X, y = phi, xout = probs, rule = 2)[[2]]
  upper = 2^approx(x = X, y = phi, xout = p.upper, rule = 2)[[2]]
  lower = 2^approx(x = X, y = phi, xout = p.lower, rule = 2)[[2]]
  results = data.frame(P, estimate, lower, upper)
  colnames(results) = c("percentile", "estimate", "lower", "upper")

  return(results)
}

#' Generate a polygon representing the confidence bounds for a grain size distribution
#'
#' \code{PolyCI} is a function that generates a data frame containing the
#' coordinates defining a polygon representing the confidence bounds around
#' a bed surface grain size distribution. The function contains an option to
#' generate a plot of the distribution showing the estimates of the percentiles
#' between the D5 and the D95, as well as the confidence limits about those
#' estimates.
#'
#' The function returns a polygon object containing y coordinates that range
#' from 0 to 1 (i.e. the proportion finer), and corresponding grain sizes that
#' define the upper and lower bounds to the grain size confidence interval.
#'
#' @param cfd is a data frame providing a list grain sizes in the first variable
#' and the corresponding cumulative proportion finer in the second. The
#' grain sizes should be recorded in mm, and the proportion finer in [0,1].
#' @param n is the total number of observations upon which the cumulative
#' frequency distribution in \code{cfd} is based
#' @param P numeric vector of percentiles defining the polygon vertices,
#'  with values in [0,100].
#' @param equaltail is a logical variable that determines whether the calculations
#' of the confidence interval are based on an approximation with equal areas in
#' each tail (the default), or based on the exact binomial solution with a coverage
#' of at least 95\%
#' @param alpha the desired confidence level for which to calculate a
#' confidence interval in [0,1].
#' @param plot optional flag to produce a graph of the resulting grain size
#' distribution
#' @export

PolyCI = function(cfd, n, P = seq(1, 99, 1), equaltail = T, alpha = 0.05, plot = FALSE){
  probs = P/100
  # use the binomial approach
  p.upper = vector(mode="numeric", length = length(probs))
  p.lower = vector(mode="numeric", length = length(probs))

  if (equaltail == T){
    for(i in seq_along(probs)){
      tmp = QuantBD(n, probs[i], alpha)
      p.upper[i] = tmp$equaltail[2]/n
      p.lower[i] = tmp$equaltail[1]/n
    }
  }else{
    for(i in seq_along(probs)){
      tmp = QuantBD(n, probs[i], alpha)
      p.upper[i] = tmp$interval[2]/n
      p.lower[i] = tmp$interval[1]/n
    }
  }

  # estimate percentiles
  phi = log2(cfd[[1]])
  X = cfd[[2]]
  estimate = 2^approx(x = X, y = phi, xout = probs, rule = 2)[[2]]
  upper = 2^approx(x = X, y = phi, xout = p.upper, rule = 2)[[2]]
  lower = 2^approx(x = X, y = phi, xout = p.lower, rule = 2)[[2]]

  x.poly = c(upper, rev(lower))
  y.poly = c(probs, rev(probs))
  poly.out = data.frame(x.poly, y.poly)

  if(plot == TRUE){
    plot(cfd[[1]], cfd[[2]],
         type = "b",
         pch = 20,
         col = "blue",
         log = "x",
         xlim = c(min(cfd[[1]]),max(cfd[[1]])),
         ylim = c(0.05,0.95),
         xlab = "grain size (mm)",
         ylab = "cum. prop. finer")
    polygon(poly.out,
            col=rgb(0, 0, 1,0.3),
            lty = 0)
    abline(h = c(0.05, 0.95))
  }
  return(poly.out)
}

#' Compare two grain size distribution samples to see if they are different
#'
#' \code{CompareCFDs} is a function that takes two cumulative frequency dist-
#' ributions with the same format as the data frames produced by MakeCFD, and
#' uses an inverse transform approach to numerically determine whether per-
#' centile estimates from the two samples are statistically different. The
#' analysis requires only the cumulative frequency distribution for each sample
#' and the number of stones measured to generate the distribution.
#'
#' The function returns a data frame listing the percentile being compared, the
#' estimate of that percentile for sample 1, the estimate for sample 2, and a
#' logical variable that indicates whether the percentiles are statistically
#' different or not.
#'
#' @param GSD1 is a data frame containing grain sizes (in mm) in the first col-
#' umn, and the proportion of the distribution finer in the second column for
#' the first sample.
#' @param GSD2 is a data frame containing grain sizes (in mm) in the first col-
#' umn, and the proportion of the distribution finer in the second column for
#' the second sample.
#' @param n1 is  the total number of observations upon which the cumulative
#' frequency distribution in \code{GSD1} is based
#' @param n2 is the total number of observations upon which the cumulative
#' frequency distribution in \code{GSD2} is based
#' @param P numeric vector of percentiles to be compared. The default is to
#' compare the 5th to the 95th percentile, in increments of 5
#' @param alpha  the desired confidence level for which to calculate a
#' confidence interval. The default is to set alpha = 0.05, for a 95\% Confidence
#' Interval
#' @export

CompareCFDs = function(GSD1, GSD2, n1, n2, P = seq(5,95, 5), alpha = 0.05){

  nr = 10^3  #number of realizations
  Dp1 = matrix(data = NA, nrow = nr, ncol = length(P))
  Dp2 = matrix(data = NA, nrow = nr, ncol = length(P))

  for (i in 1:nr) {
    # generate uniform random numbers
    u1 = runif(n1, 0, 1)
    u2 = runif(n2, 0, 1)

    # convert from uniform to shape of sediment distribution using
    # inverse transform approach
    y1 = approx(GSD1[,2], log2(GSD1[,1]), u1)[[2]]
    y1[is.na(y1)] = min(log2(GSD1[,1]))  #replicate the FINER than lim sampling effect
    y2 = approx(GSD2[,2], log2(GSD2[,1]), u2)[[2]]
    y2[is.na(y2)] = min(log2(GSD2[,1]))  #replicate the FINER than lim sampling effect

    # generate CFDs from the resampled data
    F1s = bicalc::MakeCFD(2^y1)
    F2s = bicalc::MakeCFD(2^y2)

    # interpolate Dp using inverse-transform approach
    Dp1[i,] = 2^approx(F1s$probs, log2(F1s$size), P/100)[[2]]
    Dp2[i,] = 2^approx(F2s$probs, log2(F2s$size), P/100)[[2]]
  }
  deltaDp = Dp1 - Dp2

  # hypothesis test - two tailed
  CL = matrix(data = NA, ncol = 2, nrow = length(P))
  for(i in seq_along(P)){
    CL[i,] = quantile(deltaDp[,i], c(alpha/2, 1 - alpha/2))
  }
  dffrnt = CL[,1]*CL[,2]>0

  # report sample estimates of the percentile
  D1 = 2^approx(GSD1$probs, log2(GSD1$size), P/100)[[2]]
  D2 = 2^approx(GSD2$probs, log2(GSD2$size), P/100)[[2]]

  # create a data frame presenting the results
  results = data.frame(P,D1,D2, dffrnt)
  colnames(results) = c("Percentile", "Sample_1", "Sample_2", "Stat_Diff")
  return(results)
}

#' Compare two sets of grain size measurements to see if they are different
#'
#' \code{CompareRAWs} is a function that takes two sets of grain size observations
#' and uses a resampling approach with replacement to estimate whether or not the
#' grain sizes for a percentile of interest from the two sample sets are statistically
#' different. The analysis requires the entire set of individual b-axis measurements.
#'
#' The function returns a data frame listing the percentile being compared, the
#' estimate of that percentile for sample 1, the estimate for sample 2, and a
#' logical variable that indicates whether the percentiles are statistically
#' different or not.
#'
#' @param OBS1 is a vectory containing indivuala b axis measurments (in mm) for
#' the first sample.
#' @param OBS2 is a vectory containing indivuala b axis measurments (in mm) for
#' the second sample.
#' @param P numeric vector of percentiles to be compared. The default is to
#' compare the 5th to the 95th percentile, in increments of 5
#' @param alpha  the desired confidence level for which to calculate a
#' confidence interval. The default is to set alpha = 0.05, for a 95\% Confidence
#' Interval
#' @export

CompareRAWs = function(OBS1, OBS2, P = seq(5,95, 5), alpha = 0.05){
  n1 = length(OBS1)
  n2 = length(OBS2)

  nr = 10^3  #number of realizations
  Dp1 = matrix(data = NA, nrow = nr, ncol = length(P))
  Dp2 = matrix(data = NA, nrow = nr, ncol = length(P))
  for (i in 1:nr) {
    # resample with replacement
    y1 = sample(OBS1, n1, replace = T)
    y2 = sample(OBS2, n2, replace = T)

    # extract the percentiles from the resampled data sets
    Dp1[i,] = as.numeric(quantile(y1, probs = P/100))
    Dp2[i,] = as.numeric(quantile(y2, probs = P/100))
  }
  deltaDp = Dp1 - Dp2

  # hypothesis test - two tailed
  CL = matrix(data = NA, ncol = 2, nrow = length(P))
  for(i in seq_along(P)){
    CL[i,] = quantile(deltaDp[,i], c(alpha/2, 1 - alpha/2))
  }
  dffrnt = CL[,1]*CL[,2]>0

  # report sample estimates of the percentile
  D1 = as.numeric(quantile(OBS1, probs = P/100))
  D2 = as.numeric(quantile(OBS2, probs = P/100))

  # create a data frame presenting the results
  results = data.frame(P,D1,D2, dffrnt)
  colnames(results) = c("Percentile", "Sample_1", "Sample_2", "Stat_Diff")
  return(results)
}
