# combined.R: Script to fit a (point count)+(ARU) occupancy model
#
# Usage:
#   * Run interactively
#   * or ```   source('models/src/combined.R')   ```


# libraries ---------------------------------------------------------------
library(future)
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(promises)
library(tidyverse)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib.R"))


guilds <- read_csv(
  "models/input/bird_guilds.csv",
  col_names = c("name", "code6", "code4", "guild"),
  col_types = cols(
    name = col_character(),
    code6 = col_character(),
    code4 = col_character(),
    guild = col_character()
  )
)
speciesCodes <- guilds$code4

# parameters --------------------------------------------------------------
year <- 2021
threshold <- -3
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered)

# data --------------------------------------------------------------------
# drive_auth(email = TRUE) # do not prompt when only one email has token
# drive_sync(
#  here("acoustic/data_ingest/output/"),
#  "https://drive.google.com/drive/folders/1eOrXsDmiIW9YqJWrlUWR9-Cgc7hHKD_5"
# )


runTrial <- function(speciesCode) {
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
    p10 ~ dbeta(2, 10) # p10 = Pr(y = 1 | z = 0)
    p11 ~ dbeta(5, 2) # p11 = Pr(y = 1 | z = 1)
    lam ~ dunif(0, 1000) # lambda: rate of target-species calls detected
    ome ~ dunif(0, 1000) # omega: rate of non-target detections

    # Parameters of the observation model for the scores
    mu[1] ~ dnorm(1, 1)
    mu[2] ~ dnorm(-2, 0.5)
    sigma[1] ~ dunif(0, 10)
    tau[1] <- 1 / (sigma[1] * sigma[1])
    sigma[2] ~ dunif(0, 10)
    tau[2] <- 1 / (sigma[2] * sigma[2])

    # Likelihood part 1: detection data and ARU counts
    for (i in 1:nsites) { # Loop over sites
      z[i] ~ dbern(psi) # Latent occupancy states

      p[i] <- z[i]*p11 + (1-z[i])*p10 # Detection probability including over-detections

      for(j in 1:nsurveys.pc) { # Loop over occasions
        y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
      }

      for(j in 1:nsurveys.aru) { # Loop over occasions
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
      N1[k] <- ifelse(g[k]==1, 1, 0)
    }

    # Derived quantities
    Npos <- sum(N1[])
  }
  ")

  # initialization
  zst <- rep(1, data$nsites)
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > threshold] <- 1
  gst[data$score <= threshold] <- 2
  inits <- function() {
    list(
      mu = c(rnorm(1, 1, 1), rnorm(1, -2, 0.5)),
      sigma = c(runif(1, 0, 10), runif(1, 0, 10)), z = zst,
      psi = runif(1, 0, 1), p10 = rbeta(1, 2, 10), p11 = rbeta(1, 5, 2),
      lam = runif(1, 0, 1000), ome = runif(1, 0, 1000), g = gst
    )
  }

  monitored <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "Npos")

  # MCMC settings
  na <- 2000
  ni <- 2000
  nt <- 1
  nb <- 2000
  nc <- 3

  # TODO(matt.har.vey): JAGS does not like indices in the list, since it's
  # non-numeric. Discuss team preferences on whether to omit it, nest the return
  # value of readCombined, or some other alternative.
  jagsData <- within(data, rm(indices))

  jagsResult <- jags(jagsData, inits, monitored, modelFile,
    n.adapt = na,
    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
  )
  kv <- list(result=jagsResult)
  names(kv) <- speciesCode
  return(kv)
}

trialResults <- list()
plan(multisession, workers=16)
for (speciesCode in speciesCodes) {
  promise <- future_promise({ runTrial(speciesCode) }, seed=123)
  then(
    promise,
    onFulfilled=function (kv) {
      trialResults <<- c(trialResults, kv)
    },
    onRejected=function(err) { warning(speciesCode, " failed: ", err) }
  )
}

flatMeanAndRhat <- function(jagsResult) {
  means <- unlist(jagsResult$mean)
  rhats <- unlist(jagsResult$Rhat)
  names(rhats) <- as.character(
    map(names(rhats), function(n) {
      paste("rhat", n, sep = "_")
    })
  )
  append(means, rhats)
}

getCurrentResults <- function() {
  table <- t(sapply(trialResults, flatMeanAndRhat))
  data.frame(table) %>% rownames_to_column("species")
}
