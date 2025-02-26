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

source("src/gt4ireval.R")

measure.Phi.estimate <- function(X, n_t_)
{
  g <- g.study(data = X, drop = 0)
  d <- d.study(gdata = g, queries = n_t_)
  return(d$Phi)
}

measure.Phi.actual <- function(X_estimate, X_truth)
{
  # truth
  mu_s <- colMeans(X_truth)
  mu <- mean(mu_s)
  # estimate
  mu_s.hat <- colMeans(X_estimate)

  num <- mean((mu_s - mu)^2)
  err <- mean((mu_s.hat - mu_s)^2)
  return(num / (num + err))
}
