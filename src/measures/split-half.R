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

## Compute FUN between random half splits.
## A maximum of maxTrials witll be run for each point, totaling to maxObservations at most.
split.half <- function(X, FUN, points = 20, maxObservations = 1000, maxTrials = 100)
{
  maxx <- nrow(X)
  # Compute sizes of topic subsets. By default, run at 20 points between 2 and n_t/2
  at <- unique(round(seq(1, floor(maxx / 2), length.out = points +1)))[-1]
  # If that results in only one point, also compute at subset size of 1 (that should result in sizes 1 and 2)
  if(length(at) == 1)
    at <- c(1, at)
  # By default, run at most 1000 times, and at most 100 trials per subset size.
  trials <- min(maxTrials, floor(maxObservations / length(at)))

  # Split-half
  observations <- data.frame(x = rep(at, each = trials), y = 0)
  for(i in seq_along(at)) {
    y <- replicate(trials, {
      x <- at[i] # subset size
      x1 <- sample(1:maxx, x) # topic ID for subset 1
      x2 <- sample((1:maxx)[-x1], x) # topic ID for subset 2
      res <- FUN(X[x1,, drop = FALSE], X[x2,, drop = FALSE]) # drop=FALSE to avoid turning matrix into vector for subset size of 1
      return(res)
    })
    observations$y[(i-1)*trials + 1:trials] <- y
  }
  return(observations)
}

## Fit the model y = a*x^b and return extrapolated values at newx points
extrapolate.exp1 <- function(x, y, newx)
{
  m <- lm(ly ~ lx, data = list(ly = log(y), lx = log(x)), subset = !is.infinite(ly))
  p <- predict(m, newdata = list(lx = log(newx)))

  return(pmax(0, pmin(1, exp(p))))
}
## Fit the model y = a*exp(b*x) and return extrapolated values at newx points
extrapolate.exp2 <- function(x, y, newx)
{
  m <- lm(ly ~ x, data = list(ly = log(y), x = x), subset = !is.infinite(ly))
  p <- predict(m, newdata = list(x = newx))

  return(pmax(0, pmin(1, exp(p))))
}
## Fit the model logit(y) = a*log(x)+b and return extrapolated values at newx points
extrapolate.logit <- function(x, y, newx)
{
  suppressWarnings(m <- glm(y ~ lx, data = list(y = y, lx = log(x)), family = binomial(link = "logit")))
  p <- predict(m, newdata = list(lx = log(newx)), type = "response")

  return(p)
}
