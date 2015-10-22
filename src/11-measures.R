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
args <- commandArgs(trailingOnly = TRUE)
measure.name <- args[1]
collection.name <- args[2]
assumptions <- as.list(as.logical(args[3:6]))
names(assumptions) <- c("normal", "homoscedastic", "uncorrelated", "random")

measure.path <- file.path("src/measures", paste0("measure.",measure.name,".R"))
original.path <- file.path("data", paste0(collection.name,".csv"))

if(!file.exists(measure.path) | !file.exists(original.path) | any(is.na(assumptions))) {
  stop("usage: Rscript src/11-measures.R <measure> <collection> <normal> <homoscedastic> <uncorrelated> <random>\n")
}

source("src/common.R")
source(measure.path)

# Prepare output directory
path <- file.path("scratch/11-measures", measure.name, collection.name, assumptionsToPath(assumptions))
dir.create(path, recursive = TRUE, showWarnings = FALSE)
path.simulation <- file.path("scratch/01-simulation", collection.name, assumptionsToPath(assumptions))

# EXECUTION ############################################################################################################

measure.estimate.FUN <- get(paste0("measure.", measure.name, ".estimate"))
measure.actual.FUN <- get(paste0("measure.", measure.name, ".actual"))

# Read original matrix of observations, to compute actual measure values
X <- as.matrix(read.csv(original.path))
X <- dropWorstSystems(X)

# Read simulations of different topic set sizes
for(n_t_ in .N_t) {
  # Output basic info and initialize directory for n_t_
  cat(paste(path, ": n_t' =", n_t_, "\n"))
  flush.console()

  # Prepare table of results: cols = trial, actual, n_t_=n_t,  n_t_=5, nt_10, ....
  results <- matrix(NA, ncol = 3 + length(.N_t), nrow = .N_TRIALS,
                    dimnames = list(NULL,c("trial","actual","n_t_n_t", paste0("n_t_",.N_t))))

  pb <- txtProgressBar(min = 0, max = .N_TRIALS, style = 3)
  setTxtProgressBar(pb, 0)
  for(trial in 1:.N_TRIALS) {
    set.seed(trial * n_t_) # All simulations for each trial-n_t_ have the same seed across assumptions

    Y <- as.matrix(read.csv(file.path(path.simulation, n_t_, paste0(trial, ".csv"))))
    actual <- measure.actual.FUN(Y, X)
    estimates <- measure.estimate.FUN(Y, .N_t)

    results[trial,] <- c(trial, actual, estimates[.N_t==n_t_], estimates)

    setTxtProgressBar(pb, trial)
  } # trials

  write.csv(file = file.path(path, paste0(n_t_,".csv")), row.names = FALSE, round(results, 6))
  close(pb)
} # n_t_
