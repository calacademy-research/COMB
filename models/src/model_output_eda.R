## Script to generate graphs from a `trialResults` output
## `trialResults` must be loaded into the environment
library(tidyverse)

# For a given species (one trialResults file), graph how the precision (1/var)
# of PropOcc changes as we add more data (point counts and ARUs)

species <- "MOUQ" # <------------------- Manually change the species here

parameter <- 'mean_psi'

precision <- function(values) {
  1 / var(values)
}

prec <- map(names(trialResults), 
                ~ trialResults[[.x]]$sims.list[[parameter]] %>% 
                  precision()
                  ) %>% unlist %>% 
  tibble(precision = .) %>% 
  mutate(nPC_nARU = names(trialResults)) %>% 
  separate(nPC_nARU, c("nPC", "nARU"), sep = "_") %>% 
  mutate(nPC = as.numeric(nPC), nARU = as.numeric(nARU))

# Plotting

ggplot(PropOcc_prec, aes(nARU, precision, color = as.factor(nPC))) +
  geom_point() +
  geom_smooth() + 
  labs(title = paste("Precision (1/var) of", parameter, "for", species))





