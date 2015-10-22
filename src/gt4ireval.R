## Copyright (C) 2013-2015  Julián Urbano <urbano.julian@gmail.com>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see http://www.gnu.org/licenses/.

# Taken from https://github.com/julian-urbano/GT4IREval/blob/d24546d90c9b4ccce2421350f2d6619d529cad70/gt4ireval.R
# For efficiency, we commented-out the interval estimates in the D-study.
# Also avoid warning messages unless required.

g.study <- function(data, drop = 0, showWarnings = FALSE) {
  if(!is.data.frame(data) & !is.matrix(data))
    stop("Data must be a data.frame or matrix (columns for systems and rows for queries)")

	# drop bottom results as indicated
	t <- data[,colMeans(data) >= quantile(colMeans(data), drop)]

	# sizes and means
	n.s <- ncol(t)
	means.s <- colMeans(t)
	n.q <- nrow(t)
	means.q <- rowMeans(t)
	means.all <- mean(means.s)

	# sums of squares and mean squares
	ss.s <- sum((means.s - means.all)^2) * n.q
	ss.q <- sum((means.q - means.all)^2) * n.s
	ss.e <- sum((t - means.all)^2) - ss.s - ss.q
	em.s <- ss.s / (n.s - 1)
	em.q <- ss.q / (n.q - 1)
	em.e <- ss.e / (n.s - 1) / (n.q - 1)

  # variance components
  var.s <- (em.s - em.e) / n.q
  var.q <- (em.q - em.e) / n.s
  var.e <- em.e
  if(var.s < 0){
    var.s <- 0
    if(showWarnings) {
      warning("Estimate of system variance is negative")
    }
  }
  if(var.q < 0){
    var.q <- 0
    if(showWarnings) {
      warning("Estimate of query variance is negative")
    }
  }
  if(var.e < 0){
    var.e <- 0
    if(showWarnings) {
      warning("Estimate of interaction variance is negative")
    }
  }

	res <- list(
		n.s = n.s, n.q = n.q,
		var.s = var.s, var.q = var.q, var.e = var.e,
		em.s = em.s, em.q = em.q, em.e = em.e
	)
	class(res) <- "gstudy"
	return(res)
}

d.study <- function(gdata, queries = gdata$n.q, stability = 0.95, alpha = 0.025) {
	if(length(alpha) > 1 & length(queries) > 1)
		stop("Only one of 'alpha' and 'queries' may have multiple values")
	if(length(alpha) > 1 & length(stability) > 1)
		stop("Only one of 'alpha' and 'stability' may have multiple values")

	# Point estimates for the indicated n.q
	Erho2_ <- gdata$var.s / (gdata$var.s + gdata$var.e / queries)
	n.q_Erho2 <- ceiling((stability * gdata$var.e) / (gdata$var.s * (1 - stability)))

	Phi_ <- gdata$var.s / (gdata$var.s + (gdata$var.q + gdata$var.e) / queries)
	n.q_Phi <- ceiling((stability * (gdata$var.q + gdata$var.e)) / (gdata$var.s * (1 - stability)))

#	# Interval estimates
#	df.s <- gdata$n.s - 1
#	df.q <- gdata$n.q - 1
#	df.e <- df.s * df.q
#
#	Erho2.lwr <- gdata$em.s / (gdata$em.e * qf(1-alpha, df.s, df.e))
#	Erho2.lwr <- (Erho2.lwr - 1) / gdata$n.q
#	n.q_Erho2.lwr <- ceiling(stability / (Erho2.lwr * (1 - stability)))
#	Erho2.lwr <- queries * Erho2.lwr / (1 + queries * Erho2.lwr)
#
#	Erho2.upr <- gdata$em.s / (gdata$em.e * qf(alpha, df.s, df.e))
#	Erho2.upr <- (Erho2.upr - 1) / gdata$n.q
#	n.q_Erho2.upr <- ceiling(stability / (Erho2.upr * (1 - stability)))
#	Erho2.upr <- queries * Erho2.upr / (1 + queries * Erho2.upr)
#
#	Phi.lwr <- gdata$em.s^2 - qf(1-alpha, df.s, Inf) * gdata$em.s * gdata$em.e
#	Phi.lwr <- Phi.lwr + (qf(1-alpha, df.s, Inf) - qf(1-alpha, df.s, df.e)) *
#		qf(1-alpha, df.s, df.e) * gdata$em.e^2
#	Phi.lwr <- Phi.lwr / ( (gdata$n.s - 1) * qf(1-alpha, df.s, Inf,) * gdata$em.s * gdata$em.e +
#		qf(1-alpha, df.s, df.q) * gdata$em.s * gdata$em.q )
#	Phi.lwr <- gdata$n.s * Phi.lwr / (gdata$n.s * Phi.lwr + gdata$n.q)
#	n.q_Phi.lwr <- ceiling(stability * (1 - Phi.lwr) / (Phi.lwr * (1 - stability)))
#	Phi.lwr <- queries * Phi.lwr / (1 + (queries - 1) * Phi.lwr)
#
#	Phi.upr <- gdata$em.s^2 - qf(alpha, df.s, Inf) * gdata$em.s * gdata$em.e
#	Phi.upr <- Phi.upr + (qf(alpha, df.s, Inf) - qf(alpha, df.s, df.e)) *
#		qf(alpha, df.s, df.e) * gdata$em.e^2
#	Phi.upr <- Phi.upr / ( (gdata$n.s - 1) * qf(alpha, df.s, Inf,) * gdata$em.s * gdata$em.e +
#		qf(alpha, df.s, df.q) * gdata$em.s * gdata$em.q )
#	Phi.upr <- gdata$n.s * Phi.upr / (gdata$n.s * Phi.upr + gdata$n.q)
#	n.q_Phi.upr <- ceiling(stability * (1 - Phi.upr) / (Phi.upr * (1 - stability)))
#	Phi.upr <- queries * Phi.upr / (1 + (queries - 1) * Phi.upr)

	res <- list(
		Erho2 = Erho2_, Phi = Phi_,
		n.q_Erho2 = n.q_Erho2, n.q_Phi = n.q_Phi,
#		Erho2.lwr = Erho2.lwr, Erho2.upr = Erho2.upr,
#		Phi.lwr = Phi.lwr, Phi.upr = Phi.upr,
#		n.q_Erho2.lwr = n.q_Erho2.lwr, n.q_Erho2.upr = n.q_Erho2.upr,
#		n.q_Phi.lwr = n.q_Phi.lwr, n.q_Phi.upr = n.q_Phi.upr,
		call = list(gstudy = gdata,
			queries = queries, stability = stability, alpha = alpha
		)
	)
	class(res) <- "dstudy"
	return(res)
}
