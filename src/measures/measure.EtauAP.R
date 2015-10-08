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

source("src/measures/tauAP.R")

measure.EtauAP.estimate <- function(X, n_t_)
{
  # Sort by system effectiveness so we don't have to worry about the sign
  X <- X[,order(colMeans(X), decreasing = T)]
  n_s <- ncol(X) # number of systems

  # Accumulate probabilities of no swap throughout the ranking
  S.e <- rep(0, length(n_t_))
  S.var <- S.e
  for(i in 2:n_s) {
    for(j in 1:(i-1)) {
      d <- X[,j] - X[,i]
      s <- pnorm(-mean(d) / sqrt(var(d) / n_t_))
      S.e <- S.e + (1-s) / (i-1)
      S.var <- S.var + s*(1-s) / (i-1)^2
    }
  }
  # Compute mean probability of no swap
  S.e <- S.e / (n_s-1)
  S.var <- S.var / (n_s-1)^2
  # And expected EtauAP
  EtauAP.e <- 2*S.e -1
  EtauAP.var <- 4*S.var

  return(EtauAP.e)
}

measure.EtauAP.actual <- tauAP.actual
