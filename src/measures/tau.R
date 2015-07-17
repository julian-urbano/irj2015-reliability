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

tau.estimate <- function(X, n_t_, extrapolator)
{
  observations <- split.half(X, tau.actual)
  observations <- observations[complete.cases(observations),]
  # Fit model to observations (upside-down and in [0,1])
  tau.hat <- extrapolator(x = observations$x, y = 1- (observations$y+1)/2, newx = n_t_)
  # Since we fitted upside down, reverse and go back to [-1,1]
  tau.hat <- 1- 2*tau.hat
  return(tau.hat)
}

tau.actual <- function(X_estimate, X_truth)
{
  return(cor(colMeans(X_estimate), colMeans(X_truth), method = "k"))
}