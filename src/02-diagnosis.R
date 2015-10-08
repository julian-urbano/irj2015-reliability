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

source("config/params.R")

# CONFIG ###############################################################################################################

# Check command-line parameters
args <- commandArgs(trailingOnly = T)
collection.name <- args[1]
assumptions <- as.list(as.logical(args[2:5]))
names(assumptions) <- c("normal", "homoscedastic", "uncorrelated", "random")

original.path <- file.path("data", paste0(collection.name,".csv"))

if(!file.exists(original.path) | any(is.na(assumptions))) {
  stop("usage: Rscript src/02-diagnosis.R <collection> <normal> <homoscedastic> <uncorrelated> <random>\n")
}

source("src/common.R")
source("src/simulation.R")

# Prepare output directory
path <- file.path("scratch/02-diagnosis", collection.name, assumptionsToPath(assumptions))
dir.create(path, recursive = T, showWarnings = F)
path.simulation <- file.path("scratch/01-simulation", collection.name, assumptionsToPath(assumptions))

# EXECUTION ############################################################################################################

# Read original matrix of observations, to compute actual indicator values
X <- as.matrix(read.csv(original.path))
X <- dropWorstSystems(X)
decX <- decompose.effects(X)
gX <- g.study(X)
# F_tX <- ecdf(decX$NU_t)
# F_EX <- apply(decX$E, 2, ecdf)
F_tX <- estimate.cdf(decX$NU_t, xmin = -1, xmax = 1)
F_EX <- apply(decX$E, 2, function(e) estimate.cdf(e, xmin = min(-1, min(e)), xmax = max(1, max(e))))

# Read simulations of different topic set sizes
for(n_t_ in .N_t) {
  # Output basic info and initialize directory for n_t_
  cat(paste(path, ": n_t' =", n_t_, "\n"))
  flush.console()

  # Prepare table of results. Diagnose X with itself to get number of columns
  diagnosis <- diagnose.simulation(X, X, decX, F_tX, F_EX, gX)
  results <- matrix(NA, ncol = 1+length(diagnosis), nrow = .N_TRIALS, dimnames = list(NULL,c("trial",names(diagnosis))))

  pb <- txtProgressBar(min = 0, max = .N_TRIALS, style = 3)
  setTxtProgressBar(pb, 0)
  for(trial in 1:.N_TRIALS) {
    set.seed(trial * n_t_) # All simulations for each trial-n_t_ have the same seed across assumptions, if needed

    Y <- as.matrix(read.csv(file.path(path.simulation, n_t_, paste0(trial, ".csv"))))
    diagnosis <- diagnose.simulation(X, Y, decX, F_tX, F_EX, gX)

    results[trial,] <- c(trial, diagnosis)

    setTxtProgressBar(pb, trial)
  } # trials

  write.csv(file = file.path(path, paste0(n_t_,".csv")), row.names = F, round(results, 6))
  close(pb)
} # n_t_

