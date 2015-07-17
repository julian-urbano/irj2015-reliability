# Copyright (C) 2015  Julián Urbano <urbano.julian@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

source("src/measures/split-half.R")

sensitivity_abs.estimate <- function(X, n_t_, extrapolator)
{
  observations <- split.half(X, sensitivity_abs.actual)
  observations <- observations[complete.cases(observations),]
  # Fit model to observations
  sensitivity.hat <- extrapolator(x = observations$x, y = observations$y, newx = n_t_)
  return(sensitivity.hat)
}

sensitivity_abs.actual <- function(X_estimate, X_truth, maxErrorRate = .05)
{
  means_estimate <- colMeans(X_estimate)
  means_truth <- colMeans(X_truth)

  # Compute deltas between all systems in estimate and truth
  d_estimate <- rep(NA, length(means_estimate)*(length(means_estimate)-1)/2)
  d_truth <- d_estimate
  n <- 0
  for(i in 1:(ncol(X_estimate)-1)) {
    for(j in (i+1):ncol(X_estimate)) {
      n <- n+1
      d_estimate[n] <- means_estimate[i] - means_estimate[j]
      d_truth[n] <- means_truth[i] - means_truth[j]
    }
  }
  # Check signs to see which are swaped (-1) and which are not (1)
  s <- sign(d_estimate)*sign(d_truth)
  # and sort abs(d_estimate) and s by absolute value of d_estimate
  d_estimate <- abs(d_estimate)
  o <- order(d_estimate, decreasing = T)
  d_estimate <- d_estimate[o]
  s <- s[o]

  correct <- s
  correct[correct != 1] <- 0
  correct <- cumsum(correct)
  incorrect <- s
  incorrect[incorrect != -1] <- 0
  incorrect <- cumsum(abs(incorrect))

  error.rates <- incorrect / (incorrect + correct)
  i <- which(error.rates <= maxErrorRate)
  if(length(i) == 0) {
    return(NA)
  } else {
    return(min(1, d_estimate[i[length(i)]])) # max 1 because in normal assumption it can be even larger
  }
}
