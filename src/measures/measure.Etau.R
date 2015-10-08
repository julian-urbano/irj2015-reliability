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

source("src/measures/tau.R")

measure.Etau.estimate <- function(X, n_t_)
{
  # Sort by system effectiveness so we don't have to worry about the sign
  X <- X[,order(colMeans(X), decreasing = T)]
  # Number of systems and system pairs
  n_s <- ncol(X)
  n_s.pairs <- n_s * (n_s - 1) / 2

  # Accumulate probabilities of swap through all pairs of systems
  S.e <- rep(0, length(n_t_))
  S.var <- S.e
  for(i in 1:(n_s - 1)) {
    for(j in (i + 1):n_s) {
      d <- X[,i] - X[,j]
      s <- pnorm(-mean(d) / sqrt(var(d) / n_t_))
      S.e <- S.e + s
      S.var <- S.var + s*(1-s)
    }
  }
  # Compute mean probability of swap
  S.e <- S.e / n_s.pairs
  S.var <- S.var / n_s.pairs^2
  # And expected Etau
  Etau.e <- 1 - 2*S.e
  Etau.var <- 4*S.var

  return(Etau.e)
}

measure.Etau.actual <- tau.actual
