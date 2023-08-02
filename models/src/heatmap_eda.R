## Script to make heatmap of outputs from model
library(lattice)
library(tidyverse)

# Aggregating the results
x <- getCurrentResults()

results <- expand_grid(1:3, 1:24) %>% 
  rename(nPC = 1, nARU = 2) %>% 
  left_join(x, by = c("nPC", "nARU"))

# Making a matrix of the psi for heatmapping
psi_matrix <- matrix(results$psi, nrow = 3, ncol = 24, byrow = T, 
                         dimnames = list(1:3, 1:24))

levelplot(t(psi_matrix), col.regions=heat.colors(100, rev = T), 
          xlab = "nARU", ylab = "nPC", main = list("psi parameter"))


# Making a matrix of the 95% confidence intervals for psi
psi_CI_matrix <- matrix(results$psi_0.95_CI, nrow = 3, ncol = 24, byrow = T, 
                     dimnames = list(1:3, 1:24))

levelplot(t(psi_CI_matrix), col.regions=heat.colors(100), 
          xlab = "nARU", ylab = "nPC", main = list("psi 95% Confidence Interval"))


ggplot(results, aes(x = nARU, y = psi_0.95_CI)) +
  geom_smooth() +
  geom_point()
