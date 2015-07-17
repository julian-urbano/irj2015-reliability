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

tauAP.estimate <- function(X, n_t_, extrapolator)
{
  observations <- split.half(X, tauAP.actual)
  observations <- observations[complete.cases(observations),]
  # Fit model to observations (upside-down and in [0,1])
  tauAP.hat <- extrapolator(x = observations$x, y = 1- (observations$y+1)/2, newx = n_t_)
  # Since we fitted upside down, reverse and go back to [-1,1]
  tauAP.hat <- 1- 2*tauAP.hat
  return(tauAP.hat)
}

tauAP.actual <- function(X_estimate, X_truth) {
  n_s <- ncol(X_truth) # number of systems

  truth <- colMeans(X_truth)
  estimate <- colMeans(X_estimate)
  truth.sorted <- order(truth, decreasing = T)
  estimate.sorted <- order(estimate, decreasing = T)

  p <- 0
  for(i in 2:n_s) {
    id_i <- truth.sorted[i]
    above_truth <- truth.sorted[1:i] # include i-th system, so we don't have to take care of the case i=1
    above_estimate <- estimate.sorted[1:which(estimate.sorted == id_i)]
    C_i <- sum(above_estimate %in% above_truth) -1 # -1 to remove the i-th included above
    p <- p + C_i / (i-1)
  }
  return(2*p / (n_s-1) -1)
}
