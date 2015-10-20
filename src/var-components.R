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

## Returns the effects in 'allEffects' that contain all of 'factor.names' and with a given interaction order.
## If 'factor.names' is NULL, it returns all effects of the given interaction order.
effectsOfFactors <- function(allEffects, factor.names = NULL, interaction.order = 1) {
  w <- strsplit(allEffects, ":", fixed = TRUE) # split effects in their factors
  
  w <- sapply(w, function(e) { return(
    (is.null(factor.names) | # (any factors, or 
       length(intersect(factor.names, e)) == length(factor.names)) & # include all wanted factors), and
      length(e) == interaction.order # wanted interaction order
  )})
  return(allEffects[w])
}

## Returns an ANOVA table for the model 'm.lm' extended with the estimated variance components for each effect, as well
## as the order of the interaction. If 'negativeToZero' is TRUE, negative variance estimates are replaced with 0.
varComponents <- function(m.lm, negativeToZero = TRUE) {
  # Run ANOVA and initialize var and n.factors columns
  anovaTable <- anova(m.lm)
  anovaTable$var <- NA
  anovaTable$n.factors <- 0
  
  allEffects <- rownames(anovaTable)[-nrow(anovaTable)] # names of all effects in the ANOVA table, except for residuals
  mainFactors <- effectsOfFactors(allEffects) # names of main factors
  
  nonZeroResidual <- anovaTable["Residuals", "Df"] != 0 # do we have replicates?
  
  # Iterate all effects
  for(effect in strsplit(allEffects, ":", fixed = TRUE)) {
    effect.var <- 0
    interaction.orders <- length(effect):length(mainFactors) # order of this interaction and all the subsequent ones
    
    # Iterate all orders
    for(i in seq_along(interaction.orders)) {
      interaction.effects <- effectsOfFactors(allEffects, effect, interaction.orders[i])
      if(length(interaction.effects) > 0) {
        # If odd iteration we add; if even iteration we subtract
        effect.var <- effect.var + (-1)^(i-1) *
          sum(anovaTable[interaction.effects, "Mean Sq"])
      }else{
        # If there are no interaction effects of this order, skip
        i <- i-1
        break
      }
    }
    
    # Subtract (or add) residual variance in case of replicates
    if(nonZeroResidual) {
      effect.var <- effect.var + (-1)^i * anovaTable["Residuals", "Mean Sq"]
    }
    
    effect.var <- effect.var / length(m.lm$residuals) * prod(anovaTable[effect, "Df"] +1)
    anovaTable[effectsOfFactors(allEffects, effect, length(effect)),
               c("var", "n.factors")] <- c(effect.var, length(effect))
  }
  
  # Compute residual variance in case of replicates
  if(nonZeroResidual) {
    anovaTable["Residuals", "var"] <- anovaTable["Residuals", "Mean Sq"]
  }
  
  # set negative estimates to zero if asked to
  if(negativeToZero) {
    anovaTable$var <- pmax(0, anovaTable$var)
  }
  
  # Add percentage colum
  anovaTable$var.pct <- anovaTable$var / sum(anovaTable$var)
  
  return(anovaTable)
}
