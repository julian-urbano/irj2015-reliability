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
measure.name <- args[1]

measure.path <- file.path("scratch/11-measures", measure.name)

#if(!dir.exists(measure.path)) { # requires R > 3.2.0
#  stop("usage: Rscript src/12-measures-compile.R <measure>\n")
#}

source("src/common.R")
source("src/var-components.R")

# COLLECT DATA #########################################################################################################

# Prepare output directory
path <- "output/measures"
dir.create(path, recursive = T, showWarnings = F)

data.measure <- NULL # to store measure-wise data

for(collection.name in list.dirs(measure.path, F, F)) {
  # Output basic info
  cat(paste(measure.path, ": collection =", collection.name, "\n"))
  flush.console()
  
  data.collection <- NULL # to store collection-wise data, to be appended to data.measure
  
  pb <- txtProgressBar(min = 0, max = 16 * length(.N_t), style = 3)
  setTxtProgressBar(pb, 0)
  i <- 0
  for(assumptions.name in list.dirs(file.path(measure.path, collection.name), F, F)) {
    assumptions <- strsplit(assumptions.name, split="")[[1]]
    for(n_t in .N_t) {
      data <- read.csv(file.path(measure.path, collection.name, assumptions.name, paste0(n_t, ".csv")))
      # Append info about collection, assumptions and n_t, and append to data-collection
      data <- cbind(stringsAsFactors = F,
                    collection = collection.name, normal = assumptions[1], homoscedastic = assumptions[2],
                    uncorrelated = assumptions[3], random = assumptions[4], n_t = n_t, data)
      data.collection <- rbind(data.collection, data)
      
      i <- i+1
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)
  data.measure <- rbind(data.measure, data.collection)
}

write.csv(file = file.path(path, paste0(measure.name, ".csv")), row.names = F, data.measure)

# VAR COMPONENTS #######################################################################################################

# Existing collection ==================================================================================================
cat("var components (existing)\n")
flush.console()

# Prepare output directory
path <- "output/measures/variance-existing"
dir.create(path, recursive = T, showWarnings = F)

df <- data.measure #read.csv(file.path("output/measures", paste0(measure.name, ".csv")))
df$n_t <- as.factor(df$n_t) # for ANOVA

m.lm <- lm(I(n_t_n_t-actual)~collection+n_t+random+uncorrelated+homoscedastic+normal, data = df)
m.anova <- varComponents(m.lm, negativeToZero = FALSE)
write.csv(file = file.path(path, paste0(measure.name,".csv")), m.anova)

# New collection =======================================================================================================
cat("var components (new)\n")
flush.console()

# Prepare output directory
path <- "output/measures/variance-new"
dir.create(path, recursive = T, showWarnings = F)

df <- data.measure #read.csv(file.path("output/measures", paste0(measure.name, ".csv")))

# Let's transform our data so that we end up with actual-estimated for all combinations of n_t and n_t'
# This is probably easier using some of dplyr or reshape packages...
res <- aggregate(actual~collection+normal+homoscedastic+uncorrelated+random+n_t, df, mean) # mean observed accuracies

s <- stack(df[-1:-9]) # transform into long format (skip collection, assumptions, etc)
names(s)[1:2] <- c("estimate","n_tp")
s$n_tp <- as.integer(gsub("n_t_","",s$n_tp,fixed = T)) # transform n_t' values into integers

# append all our factors; they are in the same order as in s
s$collection <- df$collection
s$normal <- df$normal
s$homoscedastic <- df$homoscedastic
s$uncorrelated <- df$uncorrelated
s$random <- df$random
s$n_t <- df$n_t
s$trial <- 1:.N_TRIALS

# sort and append observed scores
s <- s[order(s$collection, s$normal, s$homoscedastic, s$uncorrelated, s$random, s$n_tp, s$n_t, s$trial),]
res <- res[order(res$collection, res$normal, res$homoscedastic, res$uncorrelated, res$random, res$n_t),]
s$actual <- rep(res$actual, each = length(.N_t)*.N_TRIALS)

s$n_t <- as.factor(s$n_t) # for ANOVA
s$n_tp <- as.factor(s$n_tp) # for ANOVA

# And finally...fit model and decompose variance
m.lm <- lm(I(estimate-actual)~collection+n_t*n_tp+random+uncorrelated+homoscedastic+normal, data = s)
m.anova <- varComponents(m.lm, negativeToZero = FALSE)
write.csv(file = file.path(path, paste0(measure.name,".csv")), m.anova)
