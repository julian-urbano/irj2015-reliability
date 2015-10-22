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

# Note that this script will work only with the collections from the paper.
# You may reuse this code with other collections, but will need some changes for proper display.

# Pretty names of the collections.
.COLLECTION_PRETTY_NAMES <- c("Enterprise","Genomics","Robust","Web")

source("config/params.R")
source("src/common.R")
source("src/gt4ireval.R")
source("src/simulation.R")
source("src/measures/tau.R")
source("src/measures/sensitivity_rel.R")

# PLOT DEFINITIONS #####################################################################################################

my.dev.pdf <- function(num, ratio=.8, file, ...){
  width <- my.dev.width(num)
  height <- width*ratio
  pdf(file=file, width=width, height=height)
  my.dev.par(...)
}
my.dev.win <- function(num, ratio=.8, file, ...){
  width <- my.dev.width(num)
  height <- width*ratio
  dev.new(width=width, height=height)
  my.dev.par(...)
}
my.dev.width <- function(num=1){
  ret = 4;
  if(num != 1) ret = 6 / num
  return(ret*1.8)
}
# change default plot margins
my.dev.par <- function(mar=c(3.15,3.3,2,0.4), mgp=c(1.9,.7,0), ...){
  par(mar=mar, mgp=mgp, ...)
}
my.dev.abline <- function(col="darkgrey", lwd=1, lty=2, ...){
  abline(col=col, lwd=lwd, lty=lty, ...)
}
my.dev.new <- my.dev.pdf; my.dev.off <- function(...) { off <- capture.output(dev.off(...)) }
#my.dev.new <- my.dev.win; my.dev.off <- function(){}

## Used to plot just a sample of outliers in a boxplot.
plot.outlier.sample <- function(b, maxPerGroup = 100, seed = NULL, ...)
  # b : the result of calling boxplot(..., outpch = NA).
{
  if(!is.null(seed)) {
    set.seed(seed)
  }

  i <- numeric()
  for(group in unique(b$group)) {
    i_group <- which(b$group == group)
    if(length(i_group) <= maxPerGroup)
      i <- c(i, i_group)
    else
      i <- c(i, sample(i_group, size = maxPerGroup))
  }
  points(b$group[i], b$out[i], ...)
}

# FIGURE 1: Assumptions ################################################################################################

dir.create("output/paper/fig1", showWarnings = FALSE, recursive = TRUE)

# Observed and Residuals (Normal) ======================================================================================
df <- as.matrix(read.csv("data/web2004.csv"))
df <- dropWorstSystems(df, .25) # drop bottom 25% of systems
dec <- decompose.effects(df)

my.dev.new("output/paper/fig1/ass-observed.pdf", num = 3, ratio = .7)
hist(df[,"sys59"], probability = TRUE,
     main = "Observed scores", xlab = "Reciprocal Rank")
my.dev.off()

my.dev.new("output/paper/fig1/ass-residuals.pdf", num = 3, ratio = .7)
hist(dec$E[,"sys59"], probability = TRUE,
     main = "Residuals", xlab = "Reciprocal Rank")
x <- seq(-1,1,.01)
lines(x, dnorm(x, mean = 0, sd = sd(dec$E[,"sys59"])), col = "red")
my.dev.off()

# Correlations =========================================================================================================
df <- as.matrix(read.csv("data/enterprise2006.csv"))
df <- dropWorstSystems(df, .25) # drop bottom 25% of systems
dec <- decompose.effects(df)

my.dev.new("output/paper/fig1/ass-cor1.pdf", num = 3, ratio = .8)
plot(dec$NU_t, dec$E[,"sys53"],
     main = "Correlated effects", xlab = "Topic effect", ylab = "System residual")
abline(lm(dec$E[,"sys53"] ~ dec$NU_t), col = "red")
my.dev.off()

my.dev.new("output/paper/fig1/ass-cor2.pdf", num = 3, ratio = .8)
plot(dec$E[,"sys88"], dec$E[,"sys74"],
     main = "Correlated systems", xlab = "Residual of first system", ylab = "Residual of second system")
abline(lm(dec$E[,"sys74"] ~ dec$E[,"sys88"]), col = "red")
my.dev.off()

# FIGURE 2: Split-half extrapolation ###################################################################################

dir.create("output/paper/fig2", showWarnings = FALSE, recursive = TRUE)

X <- as.matrix(read.csv("data/genomics2004.csv"))
X <- dropWorstSystems(X, .25) # drop bottom 25% of systems

# tau ==================================================================================================================
set.seed(1)

my.dev.new("output/paper/fig2/split-half1.pdf", num = 3, ratio = .8)
obs <- split.half(X, FUN = tau.actual)
plot(obs, xlim = c(1,100), ylim = c(0,1), yaxs="i", xaxs="i", las=3, col = rgb(0,0,0,.1),
     main = expression("Split-Half Extrapolation of "*tau), xlab = expression(n*"'"[t]), ylab = expression(tau))
a <- aggregate(y~x, obs, mean)
points(a, pch = 19)

newx <- seq(1, 100, 1)
tau.hat <- 1- extrapolate.exp1(x = obs$x, y = 1- (obs$y+1)/2, newx = newx)*2
lines(newx, tau.hat, col = 2)
tau.hat <- 1- extrapolate.exp2(x = obs$x, y = 1- (obs$y+1)/2, newx = newx)*2
lines(newx, tau.hat, col = 3)
tau.hat <- 1- extrapolate.logit(x = obs$x, y = 1- (obs$y+1)/2, newx = newx)*2
lines(newx, tau.hat, col = 4)
legend("bottomright", lwd = 1, col = 2:4, c("exp1", "exp2", "logit"))
my.dev.off()

# relative sensitivity =================================================================================================
set.seed(1)

my.dev.new("output/paper/fig2/split-half2.pdf", num = 3, ratio = .8)
obs <- split.half(X, FUN = sensitivity_rel.actual)
plot(obs, xlim = c(1,100), ylim=c(0,1), yaxs="i", xaxs="i", las=3, col = rgb(0,0,0,.1),
     main = expression("Split-Half Extrapolation of "*sens[rel]),
     xlab =  expression(n*"'"[t]), ylab=expression(sens[rel]))
a <- aggregate(y~x, obs, mean)
points(a, pch = 19)

newx <- seq(1, 100, 1)
rsens.hat <- extrapolate.exp1(x = obs$x, y = obs$y, newx = newx)
lines(newx, rsens.hat, col = 2)
rsens.hat <- extrapolate.exp2(x = obs$x, y = obs$y, newx = newx)
lines(newx, rsens.hat, col = 3)
rsens.hat <- extrapolate.logit(x = obs$x, y = obs$y, newx = newx)
lines(newx, rsens.hat, col = 4)
my.dev.off()

# FIGURE 3: Beta sampling ##############################################################################################

dir.create("output/paper/fig3", showWarnings = FALSE, recursive = TRUE)

Z <- seq(0, 1, .01)
my.dev.new("output/paper/fig3/beta-sampling.pdf", num = 2, ratio = .6)
plot(NA, main = "Non-random sampling of topic effects",
     xlab = expression(nu[t]*" quantile"), ylab = "Probability of being sampled",
     xlim = 0:1, ylim = c(0,.1), xaxs = "i", yaxs = "i", las = 3)
shapes <- expand.grid(c(.01, 2), c(2, 8)) # alpha and beta parameters
for(s in 1:nrow(shapes)) {
  # This is more complex than needed for the purpose of the figure, but I kept it this way because this is how it is
  # done in the simulation. Take Z (the percentiles), and use |Z|+2 quantiles. Compute Beta densities, and remove the
  # first and last values, just in case they are infinity. In the end, we end up with |Z| values again to plot.
  P <- dbeta(seq(0, 1, length.out = length(Z)+2), shapes[s,1], shapes[s,2])[-c(1,length(Z)+2)]
  lines(Z, P/sum(P), col = s)
  lines(rev(Z), P/sum(P), col = s, lty = 2)
}
abline(h = 1/length(Z), col = "darkgrey")
legend("topright", lwd = 1, lty = rep(1, nrow(shapes)+1), col = c(palette()[1:nrow(shapes)], "darkgrey"), bg = "white",
       c(apply(shapes, 1, function(s) paste(c("alpha =",", beta ="), s, collapse= " ")), "Uniform sampling"))
my.dev.off()

# TABLE 1: TREC collections ############################################################################################

dir.create("output/paper/table1", showWarnings = FALSE, recursive = TRUE)

cat(file = "output/paper/table1/table.txt", append = FALSE,
    paste(c("Track", "Measure", "$n_s$", "$(n_s)$", "$n_t$",
            "$\\sigma^2_s$", "$\\sigma^2_t$", "$\\sigma^2_{st}$"),
          collapse = " & "), "\\\\","\n")
for(i in seq_along(.COLLECTION_PRETTY_NAMES)) {
  filename <- list.files("data", "*.csv", full.names = TRUE)[i]
  track.name <- basename(tools::file_path_sans_ext(filename))
  df <- as.matrix(read.csv(filename))

  n_s <- ncol(df) # number of systems with bottom 25%
  n_t <- nrow(df)
  df <- dropWorstSystems(df, .25) # drop bottom 25% of systems
  n_s2 <- ncol(df) # number of systems without bottom 25%

  g <- g.study(df, drop = 0)
  vars <- c(g$var.s, g$var.q, g$var.e)
  vars <- round(vars/sum(vars)*100)
  cat(file = "output/paper/table1/table.txt", append = TRUE,
      paste(c(track.name, "??", n_s2, paste0("(",n_s,")"), n_t, vars), collapse = " & "), "\\\\", "\n")
}

# FIGURE 4-5: Diagnosis boxplots #######################################################################################

dir.create("output/paper/fig4-5", showWarnings = FALSE, recursive = TRUE)

min_n_t <- 15
df <- read.csv("output/diagnosis/indicators.csv")
df <- df[df$n_t >= min_n_t,] # filter out small collections

set.seed(1) # for plotting random subsets of outliers

# mu_s =================================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-mus1.pdf", num = 2, ratio = .7)
b <- boxplot(mu_s~random+collection, df, outpch = NA,
             main = expression(hat(mu)[s]-mu[s]*" by Random sampling"), ylab = "Deviation",
             names = rep(c("No","Yes"), 4), ylim=c(-.5,.5))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(mu_s~random+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

my.dev.new("output/paper/fig4-5/diagnosis-mus2.pdf", num = 2, ratio = .7)
b <- boxplot(mu_s~n_t+random, df, outpch = NA,
             main = expression(hat(mu)[s]-mu[s]*" by "*n*"'"[t]), ylab = "Deviation",
             names = rep(.N_t[.N_t >= min_n_t], 2), las = 3, ylim=c(-.5,.5))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(mu_s~n_t+random, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(6,18), c("Non-random sampling", "Random sampling"), lty = 0)
my.dev.off()

# tau_AP ===============================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-tauAP2.pdf", num = 2, ratio = .7)
b <- boxplot(tauAP~n_t+random, df, outpch = NA,
             main = expression(tau[AP]*" by "*n*"'"[t]), ylab = "Correlation", names = rep(.N_t[.N_t >= min_n_t], 2),
             las = 3, ylim=c(0,1))
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(tauAP~n_t+random, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(6,18), c("Non-random sampling", "Random sampling"), lty = 0)
my.dev.off()

# w2_t =================================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-w2t1.pdf", num = 2, ratio = .7)
b <- boxplot(w2_t~random+collection, df, outpch = NA,
             main = expression(hat(omega)[t]^2*" by Random sampling"), ylab = "Distance (log-scaled)",
             names = rep(c("No","Yes"), 4), log="y")
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(w2_t~random+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

my.dev.new("output/paper/fig4-5/diagnosis-w2t2.pdf", num = 2, ratio = .7)
b <- boxplot(w2_t~n_t+random, df, outpch = NA,
             main = expression(hat(omega)[t]^2*" by "*n*"'"[t]), ylab = "Distance (log-scaled)",
             names = rep(.N_t[.N_t >= min_n_t], 2), las = 3, log="y")
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(w2_t~n_t+random, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(6,18), c("Non-random sampling", "Random sampling"), lty = 0)
my.dev.off()

# w2_E =================================================================================================================

my.dev.new("output/paper/fig4-5/diagnosis-w2E.pdf", num = 2, ratio = .7)
b <- boxplot(w2_E~n_t+random, df, outpch = NA,
             main = expression(hat(omega)[st]^2*" by "*n*"'"[t]), ylab = "Distance (log-scaled)",
             names = rep(.N_t[.N_t >= min_n_t], 2), las = 3, log="y")
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(w2_E~n_t+random, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(6,18), c("Non-random sampling", "Random sampling"), lty = 0)
my.dev.off()

# var_s ================================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-vars.pdf", num = 2, ratio = .7)
b <- boxplot(var_s~random+collection, df, outpch = NA,
             main = expression(hat(sigma)^2*(s)-sigma^2*(s)*" (%) by Random sampling"), ylab = "Deviation",
             names = rep(c("No","Yes"), 4), ylim = c(-.11,.3))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(var_s~random+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

# var_t ================================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-vart.pdf", num = 2, ratio = .7)
b <- boxplot(var_t~random+collection, df, outpch = NA,
             main = expression(hat(sigma)^2*(t)-sigma^2*(t)*" (%) by Random sampling"), ylab = "Deviation",
             names = rep(c("No","Yes"), 4), ylim = c(-.7,.2))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(var_t~random+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

# var_st ===============================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-varst.pdf", num = 2, ratio = .7)
b <- boxplot(var_st~random+collection, df, outpch = NA,
             main = expression(hat(sigma)^2*(st)-sigma^2*(st)*" (%) by Random sampling"), ylab = "Deviation",
             names = rep(c("No","Yes"), 4), ylim = c(-.2,.7))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(var_st~random+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

# sd_var ===============================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-sdvar1.pdf", num = 2, ratio = .7)
b <- boxplot(sd_var~homoscedastic+collection, df, outpch = NA,
             main = expression(sd(bold(E[s]*sigma^2*(hat(nu)[st])))-
                                 sd(bold(E[s]*sigma^2*(nu[st])))*" by Homoscedasticity"), ylab = "Deviation",
             names = rep(c("No","Yes"), 4), ylim = c(-.025,.03))
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(sd_var~homoscedastic+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

my.dev.new("output/paper/fig4-5/diagnosis-sdvar2.pdf", num = 2, ratio = .7)
b <- boxplot(sd_var~n_t+homoscedastic, df, outpch = NA,
             main = expression(sd(bold(E[s]*sigma^2*(hat(nu)[st])))-sd(bold(E[s]*sigma^2*(nu[st])))*" by n'"[t]),
             ylab = "Deviation", names = rep(.N_t[.N_t >= min_n_t], 2), ylim = c(-.025,.03), las = 3)
abline(h = 0, lty = 2, col = "darkgrey") # line at zero
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(sd_var~n_t+homoscedastic, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(6,18), c("Heteroscedastic", "Homoscedastic"), lty = 0)
my.dev.off()

# cor ==================================================================================================================
my.dev.new("output/paper/fig4-5/diagnosis-cor.pdf", num = 2, ratio = .7)
b <- boxplot(cor~uncorrelated+collection, df, outpch = NA,
             main = expression("cor"*(hat(Sigma)*","*Sigma)*" by Uncorrelated effects"), ylab = "Correlation",
             names = rep(c("No","Yes"), 4), ylim = c(0,1))
abline(h = 1, lty = 2, col = "darkgrey") # line at one
plot.outlier.sample(b, maxPerGroup = 50) # outliers

a <- aggregate(cor~uncorrelated+collection, df, FUN = mean) # compute and plot means
points(a[,3], col = "red", pch = 4, lwd = 2)

axis(1, line = 1.4, at = c(1,3,5,7)+.5, .COLLECTION_PRETTY_NAMES, lty = 0)
my.dev.off()

# TABLE 2: Diagnosis var components ####################################################################################

dir.create("output/paper/table2", showWarnings = FALSE, recursive = TRUE)

vcs <- NULL
vcs.totals <- c()
indicator.names <- c("mu_s", "tauAP", "w2_t", "w2_E", "var_s", "var_t", "var_st", "sd_var", "cor")
effect.names <- c("normal", "homoscedastic", "uncorrelated", "random", "n_t", "collection", "Residuals")

# Iterate indicators (could be easier with a package like xtable, but anyway...)
for(f in indicator.names) {
  df <- read.csv(paste0("output/diagnosis/variance/",f,".csv"), row.names = 1)
  df$var.pct <- df$var.pct*100
  # append var components and total variance to vcs
  vcs <- cbind(vcs, df[effect.names, "var.pct"])
  vcs.totals <- c(vcs.totals, sum(df$var))
}

# Mark cells wih <1, round, add the '%' to each cell, and add total variance row
vcs[vcs < 1] <- NA
vcs <- round(vcs)
vcs[is.na(vcs)] <- "$<\\!1$"
vcs[seq_along(vcs)] <- paste0(vcs, "\\%")
colnames(vcs) <- indicator.names
rownames(vcs) <- c("Normality", "Homoscedasticity", "Uncorrelated effects", "Random sampling",
                   "$n_t'$", "Collection", "residuals")
vcs <- rbind(vcs, "Total variance" = round(vcs.totals, 6))

write.table(file = "output/paper/table2/table.txt", vcs, sep = " & ", quote = FALSE)

# FIGURE 6-8: Existing collection (bias) ###############################################################################

dir.create("output/paper/fig6-8", showWarnings = FALSE, recursive = TRUE)

measures <- c("tau", "tauAP", "sensitivity_abs", "sensitivity_rel", "Erho2", "Phi", "F1")
for(i.measure in seq_along(measures)) {
  measure.name <- measures[i.measure]

  dfs <- list() # to store all data from each estimator that goes in one plot
  if(measure.name %in% c("tau", "tauAP", "sensitivity_abs", "sensitivity_rel")) { # split-half estimators
    dfs <- c(dfs, list(read.csv(file.path("output/measures", paste0(measure.name, ".exp1.csv")))))
    dfs <- c(dfs, list(read.csv(file.path("output/measures", paste0(measure.name, ".exp2.csv")))))
    dfs <- c(dfs, list(read.csv(file.path("output/measures", paste0(measure.name, ".logit.csv")))))
  }
  if(measure.name %in% c("tau", "tauAP")) # Etau and EtauAP
    dfs <- c(dfs, list(read.csv(file.path("output/measures", paste0("E", measure.name, ".csv")))))
  if(measure.name %in% c("Erho2", "Phi", "F1"))
    dfs <- c(dfs, list(read.csv(file.path("output/measures", paste0(measure.name, ".csv")))))

  # Iterate collections, making one plot for each measure-collection pair
  collections <- unique(dfs[[1]]$collection)
  for(i.collection in seq_along(collections)) {
    collection.name <- collections[i.collection]
    # these are the indices of the data for this collection and realistic assumptions
    sset = dfs[[1]]$normal == FALSE & dfs[[1]]$homoscedastic == FALSE & dfs[[1]]$uncorrelated == FALSE & dfs[[1]]$random == TRUE &
      dfs[[1]]$collection == collection.name

    # Aggregate data from each estimator in dfs
    a <- lapply(dfs, function(df) { aggregate(I(n_t_n_t-actual)~n_t, df, mean, subset = sset) })

    # Make an empty plot for now...
    my.dev.new(file.path("output/paper/fig6-8", paste0("existing-", measure.name, "-", collection.name, ".pdf")),
               num = 4, ratio = 1)
    plot(NA, xaxs = "i", log = "x", xlim = range(.N_t), las = 3, xlab = expression(n[t]),
         main = .COLLECTION_PRETTY_NAMES[i.collection], ylim = list(c(-.25,.25),
                                                                    c(-.25,.25),
                                                                    c(-.05,.25),
                                                                    c(-.05,.25),
                                                                    c(-.25,.25),
                                                                    c(-.25,.25),
                                                                    c(-1,0))[[i.measure]],
         ylab = list(expression(hat(tau)-tau),
                     expression(hat(tau)[AP]-tau[AP]),
                     expression(widehat(sens)[abs]-sens[abs]),
                     expression(widehat(sens)[rel]-sens[rel]),
                     expression(E*hat(rho)^2-E*rho^2),
                     expression(hat(Phi)-Phi),
                     expression(hat(F)-F*" (power)"))[[i.measure]]
    )
    grid(equilogs = FALSE)
    abline(h = 0)

    # ...and here plot a line for each aggregated estimator in a
    sapply(seq_along(a), function(i) lines(a[[i]], col = i+1))

    # Finally, make the legend if this is the first collection
    if(i.collection == 1) {
      if(measure.name %in% c("tau","tauAP")){
        legend("topright", col=2:5, lwd=1, lty=1, seg.len=1.5, y.intersp = .8, bg = "white", ncol = 2,
               list(expression(exp1,exp2,logit,E*tau),
                    expression(exp1,exp2,logit,E*tau[AP]))[[i.measure]])
      }
      if(measure.name %in% c("sensitivity_abs", "sensitivity_rel")) {
        legend("topright", col=2:4, lwd=1, lty=1, seg.len=1.5, y.intersp = .8, bg = "white",
               c("exp1","exp2","logit"))
      }
    }
    my.dev.off()
  }
}

# TABLE 3-5: Existing collection (robustness) ##########################################################################

dir.create("output/paper/table3-5", showWarnings = FALSE, recursive = TRUE)

estimator.names <- list(c("tau.exp1", "tau.exp2", "tau.logit", "Etau",
                          "tauAP.exp1", "tauAP.exp2", "tauAP.logit", "EtauAP"),
                        c("sensitivity_abs.exp1", "sensitivity_abs.exp2", "sensitivity_abs.logit",
                          "sensitivity_rel.exp1", "sensitivity_rel.exp2", "sensitivity_rel.logit"),
                        c("Erho2", "Phi", "F1"))
effect.names <- c("normal", "homoscedastic", "uncorrelated", "random", "n_t", "collection", "Residuals")

# Iterate figures
for(i in seq_along(estimator.names)) {
  vcs <- NULL
  vcs.totals <- c()

  # Iterate estimators
  for(f in estimator.names[[i]]) {
    df <- read.csv(paste0("output/measures/variance-existing/",f,".csv"), row.names = 1)
    df$var.pct <- df$var.pct*100
    # append var components and total variance to vcs
    vcs <- cbind(vcs, df[effect.names, "var.pct"])
    vcs.totals <- c(vcs.totals, sum(df$var))
  }
  # Mark cells wih <1, round, add the '%' to each cell, and add total variance row
  vcs[vcs < 1] <- NA
  vcs <- round(vcs)
  vcs[is.na(vcs)] <- "$<\\!1$"
  vcs[seq_along(vcs)] <- paste0(vcs, "\\%")
  colnames(vcs) <- estimator.names[[i]]
  rownames(vcs) <- c("Normality", "Homoscedasticity", "Uncorrelated effects", "Random sampling",
                     "$n_t$", "Collection", "residuals")
  vcs <- rbind(vcs, "Total error variance" = round(vcs.totals, 6))

  write.table(file = paste0("output/paper/table3-5/table", 2+i, ".txt"), vcs, sep = " & ", quote = FALSE)
}

# FIGURE 9-11: New collection (bias) ###################################################################################

dir.create("output/paper/fig9-11", showWarnings = FALSE, recursive = TRUE)

estimator.names <- c("tau.exp1", "tau.exp2", "tau.logit", "Etau",
                     "tauAP.exp1", "tauAP.exp2", "tauAP.logit", "EtauAP",
                     "sensitivity_abs.exp1", "sensitivity_abs.exp2", "sensitivity_abs.logit",
                     "sensitivity_rel.exp1", "sensitivity_rel.exp2", "sensitivity_rel.logit",
                     "Erho2", "Phi", "F1")
for(i.estimator in seq_along(estimator.names)) {
  estimator.name <- estimator.names[i.estimator]

  df <- read.csv(file.path("output/measures", paste0(estimator.name,".csv")))

  # Here we need to work up our data a little.
  # Prepare a global matrix M to store the expected error over all collections: rows are n_t and columns are n_t'
  M <- matrix(0, nrow = length(.N_t), ncol = length(.N_t), dimnames = list(.N_t, .N_t)) # n_t x n_t'

  collections <- unique(df$collection)
  # Iterate all collections: compute expected error and add to global matrix M
  for(i.collection in seq_along(collections)) {
    collection.name <- collections[i.collection]
    # these are the indices of the data for this collection and realistic assumptions
    sset = df$normal == FALSE & df$homoscedastic == FALSE & df$uncorrelated == FALSE & df$random == TRUE &
      df$collection == collection.name

    # Let's use a temp matrix M2 for this collection: n_t x n_t' as well
    M2 <- matrix(NA, nrow = length(.N_t), ncol = length(.N_t), dimnames = list(.N_t, .N_t))
    for(i in 1:ncol(M2)) { # For each n_t'
      # get the expected accuracy for this n_t' column, given all initial n_t
      M2[,i] <- aggregate(as.formula(paste0("n_t_", .N_t[i] ,"~n_t")), df, mean, subset = sset)[,2]
    }
    # Add error to global matrix M: expected accuracy in M2 - actual accuracy
    actual <- aggregate(actual~n_t, df, mean, subset = sset)[,2]
    M <- M + (M2 - matrix(actual, ncol = ncol(M2), nrow = nrow(M2), byrow = TRUE) )
  }
  # Finally, average over collections
  M <- M / length(collections)

  my.dev.new(file.path("output/paper/fig9-11", paste0("new-", gsub(".", "-", estimator.name, fixed = TRUE),".pdf")),
             num = 4, ratio = 1)
  plot(NA, xaxs = "i", log = "x", xlim = range(.N_t), las = 3, xlab = expression(n*"'"[t]),
       main=list(expression(tau*" exp1"), expression(tau*" exp2"), expression(tau*" logit"), expression(E*tau),
                 expression(tau[AP]*" exp1"), expression(tau[AP]*" exp2"), expression(tau[AP]*" logit"), expression(E*tau[AP]),
                 expression(sens[abs]*" exp1"), expression(sens[abs]*" exp2"), expression(sens[abs]*" logit"),
                 expression(sens[rel]*" exp1"), expression(sens[rel]*" exp2"), expression(sens[rel]*" logit"),
                 expression(E*rho^2), expression(Phi), expression(F[1]))[[i.estimator]],
       ylim = list(c(-.5,.3),c(-.5,.3),c(-.5,.3),c(-.5,.3),
                   c(-.5,.3),c(-.5,.3),c(-.5,.3),c(-.5,.3),
                   c(-.2,.4),c(-.2,.4),c(-.2,.4),
                   c(-.2,.4),c(-.2,.4),c(-.2,.4),
                   c(-.2,.2),c(-.2,.2),c(-1,0))[[i.estimator]],
       ylab = expression(hat(R)-R)
  )

  grid(equilogs = FALSE)
  abline(h = 0)

  # Plot error when n_t = n_t'
  lines(.N_t, diag(M), lty = 2, lwd = 2)
  # and now the general case. Plot only a subset of initial collection sizes, so we don't clutter the plots too much
  mycols <- c(1, 1, palette()[2:6], "purple") # these are the colors we'll use to plot
  n_t <- c(5, 10, 20, 50, 100, 200)
  sapply(seq_along(n_t), function(i) { lines(.N_t, M[.N_t == n_t[i],], col = mycols[i+2]) } )

  # Finally, make the legend if needed
  if(estimator.name %in% c("tau.exp1", "tauAP.exp1", "sensitivity_abs.exp1", "sensitivity_rel.exp1", "Erho2")) {
    legend(c("topright", NA,NA,NA,
             "topright", NA,NA,NA,
             "bottomright", NA,NA,
             "bottomright", NA,NA,
             "topright", NA,NA)[i.estimator],
           lwd = c(0,1,1,0,1,1,2,1,1), lty = c(0,1,1,0,1,1,2,1,1), col = c(NA, mycols[c(3,6,1,4,7,2,5,8)]),
           ncol = 3, seg.len = 1.5, y.intersp = .8, cex = .65, x.intersp = .1, box.lty = 0, bg = "white",
           c("", expression(n[t]==5), expression(n[t]==50),
             "", expression(n[t]==10), expression(n[t]==100),
             expression(n*"'"[t]==n[t]), expression(n[t]==20), expression(n[t]==200))
    )
    box() # re-draw the box because the legend's background might hide it a little
  }

  my.dev.off()
}

# TABLE 6-9: New collection (robustness) ===============================================================================

dir.create("output/paper/table6-8", showWarnings = FALSE, recursive = TRUE)

estimator.names <- list(c("tau.exp1", "tau.exp2", "tau.logit", "Etau",
                          "tauAP.exp1", "tauAP.exp2", "tauAP.logit", "EtauAP"),
                        c("sensitivity_abs.exp1", "sensitivity_abs.exp2", "sensitivity_abs.logit",
                          "sensitivity_rel.exp1", "sensitivity_rel.exp2", "sensitivity_rel.logit"),
                        c("Erho2", "Phi", "F1"))
effect.names <- c("normal", "homoscedastic", "uncorrelated", "random",
                  "n_tp", "n_t", "n_t:n_tp", "collection", "Residuals")

# Iterate figures
for(i in seq_along(estimator.names)) {
  vcs <- NULL
  vcs.totals <- c()

  # Iterate estimators
  for(f in estimator.names[[i]]) {
    df <- read.csv(paste0("output/measures/variance-new/",f,".csv"), row.names = 1)
    df$var.pct <- df$var.pct*100
    # append var components and total variance to vcs
    vcs <- cbind(vcs, df[effect.names, "var.pct"])
    vcs.totals <- c(vcs.totals, sum(df$var))
  }
  # Mark cells wih <1, round, add the '%' to each cell, and add total variance row
  vcs[vcs < 1] <- NA
  vcs <- round(vcs)
  vcs[is.na(vcs)] <- "$<\\!1$"
  vcs[seq_along(vcs)] <- paste0(vcs, "\\%")
  colnames(vcs) <- estimator.names[[i]]
  rownames(vcs) <- c("Normality", "Homoscedasticity", "Uncorrelated effects", "Random sampling",
                     "$n_t'$", "$n_t$", "$n_t:n_t'$", "Collection", "residuals")
  vcs <- rbind(vcs, "Total error variance" = round(vcs.totals, 6))

  write.table(file = paste0("output/paper/table6-8/table", 5+i, ".txt"), vcs, sep = " & ", quote = FALSE)
}

