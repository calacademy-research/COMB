# combined.R: Script to fit a (point count)+(ARU) occupancy model
#
# Usage:
#   * Run interactively
#   * or ```   source('models/src/combined.R')   ```


# libraries ---------------------------------------------------------------
library(chron)
library(coda)
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(reshape2)
library(tidyverse)
library(wesanderson)

source(here("comb_functions.R"))


# parameters --------------------------------------------------------------
speciesCode <- "RBNU" # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- 0.5
aruVisitLimit <- 60 # only consider this many ARU visits per site (ordered)


# data --------------------------------------------------------------------
drive_auth(email = TRUE) # do not prompt when only one email has token
drive_sync(
  here("acoustic/data_ingest/output/"),
  "https://drive.google.com/drive/folders/1eOrXsDmiIW9YqJWrlUWR9-Cgc7hHKD_5"
)

# point count data frame
#
# If you want different specifications on detection distance, years included,
# etc. for your model, rerun delintPC.R with any necessary changes to those
# filters.
dataPC <- read_csv(
  here("point_counts/data_ingest/output/PC_delinted.csv"),
  col_types = cols(
    observer_fk = col_character(),
    birdCode_fk = col_character(),
    abun = col_integer(),
    point_ID_fk = col_integer(),
    year = col_integer(),
    visit = col_integer(),
    DateTime = col_datetime()
  )
)

# ARU data frame (filtered by start time)
dataML <- read_csv(
  here("acoustic/data_ingest/output/dataML_model.csv"),
  col_types = cols(
    filename = col_character(),
    File_ID = col_character(),
    ARU_ID = col_character(),
    point = col_integer(),
    Start_Time = col_datetime(),
    rebnut = col_double(),
    P = col_double(),
    Date_Time = col_datetime(),
  )
) %>%
  # morning hours only
  filter(hour(Date_Time) < 10) %>%
  # specific start minutes
  filter(minute(Date_Time) == 30 | minute(Date_Time) == 00)


# indexing ----------------------------------------------------------------
pointToIndex <- data.frame(
  point = union(
    unique(dataPC$point_ID_fk),
    unique(dataML$point)
  )
) %>%
  arrange(point) %>%
  mutate(pointIndex = seq_along(point))

visitToIndex <- dataML %>%
  select(point, Date_Time) %>%
  distinct() %>%
  group_by(point) %>%
  mutate(visit = seq_along(Date_Time))


# reshaping ---------------------------------------------------------------
nsites <- max(pointToIndex$pointIndex)
nsurveys.pc <- max(dataPC$visit)
nsurveys.aru <- aruVisitLimit

sparsePointCounts <- dataPC %>%
  filter(birdCode_fk == speciesCode, year == year) %>%
  group_by(point_ID_fk) %>%
  select(point_ID_fk, visit, abun) %>%
  inner_join(pointToIndex, by = c("point_ID_fk" = "point"))

sparseARUScores <- dataML %>%
  full_join(visitToIndex) %>%
  filter(visit <= aruVisitLimit) %>%
  inner_join(pointToIndex)

sparseARUCounts <- sparseARUScores %>%
  mutate(is.det = (rebnut > threshold)) %>%
  group_by(pointIndex, visit) %>%
  summarise(count = sum(is.det), .groups = "keep")

sparseScores <- sparseARUScores %>% filter(rebnut > threshold)
siteID <- sparseScores$pointIndex
occID <- sparseScores$visit
score <- sparseScores$rebnut
nsamples <- nrow(sparseScores)

sparseToDense <- function(rowIndices, colIndices, values, nrows, ncols) {
  m <- matrix(NA, nrows, ncols)
  for (u in 1:length(rowIndices)) {
    m[rowIndices[u], colIndices[u]] <- values[u]
  }
  m
}

# Point counts dense matrix
y.pc <- sparseToDense(
  sparsePointCounts$pointIndex,
  sparsePointCounts$visit, sparsePointCounts$abun,
  nsites,
  nsurveys.pc
)

# Point counts dense binary matrix (for Bernoulli likelihood factors)
y.ind <- (y.pc > 0) * 1

# ARU counts dense matrix
y.aru <- sparseToDense(
  sparseARUCounts$pointIndex,
  sparseARUCounts$visit,
  sparseARUCounts$count,
  nsites,
  nsurveys.aru
)


# JAGS structuring --------------------------------------------------------
data <- list(
  nsites = nsites,
  nsurveys.pc = nsurveys.pc,
  nsurveys.aru = nsurveys.aru,
  y.ind = y.ind,
  y.aru = y.aru,
  siteid = siteID,
  occid = occID,
  score = score,
  nsamples = nsamples
)


# JAGS specification ------------------------------------------------------
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
  mu[1] ~ dnorm(0, 0.01)
  mu[2] ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10)
  tau <- 1 / (sigma * sigma)

  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    z[i] ~ dbern(psi) # Latent occupancy states
    p[i] <- z[i]*p11 + (1-z[i])*p10 # Detection probability
    site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species)
    for(j in 1:nsurveys.pc) { # Loop over occasions
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    for(j in 1:nsurveys.aru) { # Loop over occasions
      y.aru[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
    }
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    # Sample specific covariate
    score[k] ~ dnorm(mu[g[k]], tau) # parameters are group specific
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
zst <- rep(1, nrow(y.pc))
gst <- sample(1:2, length(score), replace = TRUE)
gst[score > threshold] <- 1
gst[score <= threshold] <- 2
inits <- function() {
  list(
    mu = c(1, -1), sigma = 0.2, z = zst,
    psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
    lam = runif(1, 1, 2), ome = runif(1, 0, 0.4), g = gst
  )
}


# JAGS execution ----------------------------------------------------------

monitored <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "Npos")

# MCMC settings
na <- 1000
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

jagsResult <- jags(data, inits, monitored, modelFile,
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
)
