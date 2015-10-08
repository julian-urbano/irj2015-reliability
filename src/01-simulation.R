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

collection.path <- file.path("data", paste0(collection.name,".csv"))

if(!file.exists(collection.path) | any(is.na(assumptions))) {
  stop("usage: Rscript src/01-simulation.R <collection> <normal> <homoscedastic> <uncorrelated> <random>\n")
}

source("src/common.R")
source("src/simulation.R")

# Prepare output directory
path <- file.path("scratch/01-simulation", collection.name, assumptionsToPath(assumptions))
dir.create(path, recursive = T, showWarnings = F)

# EXECUTION ############################################################################################################

original <- as.matrix(read.csv(collection.path))
original <- dropWorstSystems(original)

cfg <- configure.simulation(original, assumptions$normal, assumptions$homoscedastic, assumptions$uncorrelated,
                            assumptions$random)

# Simulate different topic set sizes
for(n_t_ in .N_t) {
  # Output basic info and initialize directory for n_t_
  cat(paste(path, ": n_t' =", n_t_, "\n"))
  flush.console()
  dir.create(file.path(path, n_t_), recursive = T, showWarnings = F)

  # Run trials
  pb <- txtProgressBar(min = 0, max = .N_TRIALS, style = 3)
  setTxtProgressBar(pb, 0)
  for(trial in 1:.N_TRIALS) {
    set.seed(trial * n_t_) # All simulations for each trial-n_t_ have the same seed across assumptions
    Y <- simulate.collection(cfg, n_t_)

    write.csv(file = file.path(path, n_t_, paste0(trial, ".csv")), row.names = F, round(Y,4))
    setTxtProgressBar(pb, trial)
  } # trial
  close(pb)
} # n_t_

