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

library(pwr)
source("src/gt4ireval.R")

measure.F1.estimate <- function(X, n_t_, minDelta = .05, sig.level = 0.05)
{
  g <- g.study(X)
  sigma2 <- g$var.e
  f <- minDelta / sqrt(sigma2) * sqrt(1 / (2*ncol(X)))
  power <- pwr.f2.test(u = ncol(X)-1, v = ncol(X) * (n_t_-1), sig.level = sig.level, f2 = f^2)$power
  return(power)
}

measure.F1.actual <- function(X_estimate, X_truth, minDelta = .05, sig.level = 0.05) {
  g <- g.study(X_estimate)
  mu.hat <- colMeans(X_estimate)
  p <- pf(g$em.s / g$em.e, df1 = g$n.s-1, df2 = (g$n.s-1) * (g$n.q-1), lower.tail = FALSE)
 
  # significant and at least minDelta?
  if(p < sig.level & max(mu.hat)-min(mu.hat) >= minDelta) {
    return(1)
  } else {
    return(0)
  }
}
