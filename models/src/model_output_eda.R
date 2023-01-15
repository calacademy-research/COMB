## Script to generate graphs from a `trialResults` output
## `trialResults` must be loaded into the environment
library(tidyverse)

# For a given species (one trialResults file), graph how the precision (1/var)
# of PropOcc changes as we add more data (point counts and ARUs)

species <- "psi = 0.2, p11 = 0.4, p_aru11 = 0.4,\np_aru01 = 0.05, mu = c(-2, 1), sigma = c(0.8, 4)" # <------------------- Manually change the species here

parameter <- 'mean_psi'

precision <- function(values) {
  1 / var(values)
}

prec <- map(names(trialResults), 
                ~ trialResults[[.x]]$sims.list[[parameter]] %>% 
                  precision()
                  ) %>% unlist %>% 
  tibble(precision = .) %>% 
  mutate(params = names(trialResults)) %>% 
  separate(params, c("nPC", "nARU", "psi", "p11", "p_aru11", "p_aru01"), sep = "_") %>% 
  mutate(nPC = as.numeric(nPC), nARU = as.numeric(nARU), psi = as.numeric(psi), 
         p11 = as.numeric(p11), p_aru11 = as.numeric(p_aru11), 
         p_aru01 = as.numeric(p_aru01))

# Plotting

ggplot(prec, aes(nARU, precision, color = as.factor(nPC))) +
  geom_point() +
  geom_smooth() + 
  labs(title = paste("Precision (1/var) of", parameter, "for", species))

# Trying to find the asymptote using math

#' find_asymptote
#' @param df dataframe with precision, nPC, and nARU as cols
#' @param nPCs number of PCs to filter by in the df
#' @return the nARU where the graph asymptotes

find_asymptote <- function(df, nPCs) {
  asymptote <- NULL
  for(i in 1:24) {
    model <- lm(precision ~ nARU, data = filter(df, nPC == nPCs)[i:24,]) %>% 
      summary()
    slope <- model$coefficients[2,1] # the slope for the linear model
    if (slope < 0.25) {
      asymptote = i
      break()
    }
  }
  return(asymptote)
}
