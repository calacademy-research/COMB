# combined.R: Script to fit a (point count)+(ARU) occupancy model
#
# Usage:
#   * Run interactively
#   * or ```   source('models/src/combined.R')   ```


# libraries ---------------------------------------------------------------
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib.R"))


# parameters --------------------------------------------------------------
speciesCode <- "HAWO" # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- 0.5
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered)


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
  years = c(year),
  beginTime = dhours(6),
  endTime = dhours(10),
  visitLimit = aruVisitLimit,
  visitAggregation = "file",
  thresholdOptions = list(
    value = threshold,
    is.quantile = F
  ),
  squeeze = T
)


# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(file = modelFile, "
model {

  # Priors
  psi ~ dunif(0, 1) # psi = Pr(Occupancy)
#  p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)

  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
  }
}
")

# initialization
zst <- rep(1, data$nsites)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > threshold] <- 1
gst[data$score <= threshold] <- 2
inits <- function() {
  list(
    z = zst,
    psi = runif(1), p11 = runif(1, 0.5, 0.8)
  )
}


# JAGS execution ----------------------------------------------------------

monitored <- c("psi", "p11")

# MCMC settings
na <- 1000
ni <- 4000
nt <- 1
nb <- 1000
nc <- 6

# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data, rm(indices))

set.seed(123)

jagsResult <- jags(jagsData, inits, monitored, modelFile,
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
)
