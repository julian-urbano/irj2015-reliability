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
source("src/common.R")
source("src/var-components.R")

# COLLECT DATA #########################################################################################################

diagnosis.path <- "scratch/02-diagnosis"

# Prepare output directory
path <- "output/diagnosis"
dir.create(path, recursive = T, showWarnings = F)

data.diagnosis <- NULL # to store all data

for(collection.name in list.dirs(diagnosis.path, F, F)) {
  # Output basic info
  cat(paste(diagnosis.path, ": collection =", collection.name, "\n"))
  flush.console()

  data.collection <- NULL # to store collection-wise data, to be appended to data.diagnosis

  pb <- txtProgressBar(min = 0, max = 16*length(.N_t), style = 3)
  setTxtProgressBar(pb, 0)
  i <- 0
  for(assumptions.name in list.dirs(file.path(diagnosis.path, collection.name), F, F)) {
    assumptions <- strsplit(assumptions.name, split="")[[1]]
    for(n_t in .N_t) {
      data <- read.csv(file.path(diagnosis.path, collection.name, assumptions.name, paste0(n_t,".csv")))
      # Append info about collection, assumptions and n_t, and append to data.collection
      data <- cbind(stringsAsFactors = F,
                    collection = collection.name, normal = assumptions[1], homoscedastic = assumptions[2],
                    uncorrelated = assumptions[3], random = assumptions[4], n_t = n_t, data)
      data.collection <- rbind(data.collection, data)
      
      i <- i+1
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)
  data.diagnosis <- rbind(data.diagnosis, data.collection)
}

write.csv(file = file.path(path, "indicators.csv"), row.names = F, data.diagnosis)

# VAR COMPONENTS #######################################################################################################

# Read and filter data
df <- data.diagnosis #read.csv(file.path(path, "indicators.csv"))
df <- df[df$n_t >= 15,] # filter out small collections
df$n_t <- as.factor(df$n_t) # for ANOVA

# Prepare output directory
path <- "output/diagnosis/variance"
dir.create(path, recursive = T, showWarnings = F)

for(indicator in names(df)[-1:-7]) { # Ignore first columns (collection, assumptions, etc)
  cat(indicator,"\n")
  flush.console()

  m.lm <- lm(as.formula(paste(indicator, "~collection+n_t+random+uncorrelated+homoscedastic+normal")), data = df)
  m.anova <- varComponents(m.lm, negativeToZero = FALSE)
  write.csv(file = file.path(path, paste0(indicator,".csv")), m.anova)
}
