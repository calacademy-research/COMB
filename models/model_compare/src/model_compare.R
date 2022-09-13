# model_compare.R
# Script to change the amount of data that is ingested into model
# and compare the "confidence" and how well the model is able to 
# work with limiting the amount of data it can use

# Loading Libraries ----------------
library(tidyverse)
library(data.table)
library(here)
library(fs)
library(lubridate)
library(furrr)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib.R"))
# Ingesting data
# ATM this script should be run after the data has been ingested with combined.R
data_slice <- function(data_list, pc_days, aru_visits){
pc_days <- pc_days # number of pc days to subset
aru_visits <- aru_visits # number of ARU days to subset (possibly add amount of time per day to subset later)

pc_sub <- data_list[["indices"]][["visit.pc"]] %>% 
  filter(Visit <= pc_days)

aru_sub <- data_list[["indices"]][["visit.aru"]] %>% 
  filter(Visit_Index <= aru_visits)

data_list[["indices"]][["visit.aru"]] <- aru_sub

data_list[["indices"]][["visit.pc"]] <- pc_sub

return(data_list)
}

# Defining all possible combinations of days and visits
combinations <- expand_grid(1:3, 1:24) %>% 
  rename("pc" = 1, "aru" = 2)

all_data <- map2(.x = combinations$pc, .y = combinations$aru, 
     ~ data_slice(data, .x, .y))

model <- function(data_list){
  data <- data_list
### BEGINNING OF MODEL COPY-PASTED
modelFile <- tempfile()
cat(file = modelFile, "
model {
  # Priors
  psi ~ dunif(0, 1) # psi = Pr(Occupancy)
  p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
  p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0, 1000) # lambda: rate of target-species calls detected
  ome ~ dunif(0, 1000) # omega: rate of non-target detections
  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-1.3, 1)T(0.5,)     # -1.3, sig 1
  mu[2] ~ dnorm(-1.75, 0.1)T(0.5,)    # -1.75, sig 0.125
  sigma[1] ~ dunif(0, 10)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0, 10)
  tau[2] <- 1 / (sigma[2] * sigma[2])
  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
    # Point count
    p[i] <- z[i]*p11 + (1-z[i])*p10 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    # ARU
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
    }
  }
  # Likelihood part 2: feature score data
  for (i in 1:nsites) {
    site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species)
  }
  for(k in 1:nsamples) {
    # Sample specific covariate
    score[k] ~ dnorm(mu[g[k]], tau[g[k]]) # parameters are group specific
    probs[k,1] <- site.prob[siteid[k]]
    probs[k,2] <- 1 - site.prob[siteid[k]] # the prior class probabilities
    g[k] ~ dcat(probs[k,])
     #N1[k] <- ifelse(g[k]==1, 1, 0)
  }
  # Derived quantities
   #Npos <- sum(N1[])
}
")

zst <- rep(1, data$nsites)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > threshold] <- 1
gst[data$score <= threshold] <- 2
inits <- function() {
  list(
    mu = c(1, 0.6), sigma = c(1, 0.1), z = zst,
    psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
    lam = runif(1, 1, 2), ome = runif(1, 0, 0.4), g = gst
  )
}


# JAGS execution ----------------------------------------------------------

monitored <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "Npos")

# MCMC settings
na <- 1000
ni <- 4000
nt <- 1
nb <- 1000
nc <- 2

# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data, rm(indices))

set.seed(123)

jagsResult <- jags(jagsData, inits, monitored, modelFile,
                   n.adapt = na,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
)

return(jagsResult)
### END OF MODEL COPY-PASTED
}

plan(multisession, workers = 10)

big_jags_out <- future_map(.x = 1:72, 
                    ~ model(all_data[[.x]]))
