
# libraries ---------------------------------------------------------------
library(tidyverse)
library(chron)
library(jagsUI)
library(coda)
library(lubridate)
library(googledrive)
library(here)
library(wesanderson)
library(reshape2)


# data --------------------------------------------------------------------

drive_auth()
if (dir.exists(here("acoustic/prototype_data/")) == FALSE) {
  dir.create(here("acoustic/prototype_data/"))
}

if (dir.exists(here("acoustic/prototype_data/output/")) == FALSE) {
  dir.create(here("acoustic/prototype_data/output/"))
}

drive_sync(here("acoustic/prototype_data/output/"), "https://drive.google.com/drive/u/1/folders/1YhB2NgXqo8JLM9o7JUDqS0tHWYaBtL6S")

source(here("COMB_functions.R"))

dfc <- fread(here("point_counts/data_ingest/output/PC_delinted.csv")) # make sure these are the specifications on detection distance, years included, etc. you want for your model and rerun delintPC.R with changes to those filters if necessary

# target: y matrix [1:i, 1:j] for species HEWA and year 2019
# point count data
yRBNU <- dfc %>%
  filter(birdCode_fk == "RBNU", year == 2021) %>%
  group_by(point_ID_fk) %>%
  select(point_ID_fk, visit, abun) %>%
  inner_join(pointList, by = c("point_ID_fk" = "point"))

# ACOUSTIC DATA

dataML_model <- read_csv(here("acoustic/data_ingest/output/dataML_model.csv"))
morningML <- dataML_model %>%
  filter(hour(Date_Time) < 10) %>%
  filter(minute(Date_Time) == 30 | minute(Date_Time) == 00)

# point index
pointList <- data.frame(point = as.numeric(union(unique(dfc$point_ID_fk), unique(dataML_model$point)))) %>%
  arrange(pointList, point) %>%
  mutate(pointIndex = seq_along(point))
# filter to morning hours
visitindex <- morningML %>%
  select(point, Date_Time) %>%
  distinct() %>%
  group_by(point) %>%
  mutate(visit = seq_along(Date_Time))

visitLimit <- 60
MLscores <- morningML %>%
  full_join(visitindex) %>%
  filter(visit <= visitLimit) %>%
  inner_join(pointList)

MLcounts <- MLscores %>%
  mutate(is.det = (rebnut > 0.5)) %>%
  group_by(pointIndex, visit) %>%
  summarise(count = sum(is.det))


samples <- MLscores %>% filter(rebnut > 0.5)
siteID <- samples$pointIndex
occID <- samples$visit
nsamples <- nrow(samples)
nsites <- max(pointList$pointIndex)
nsurveys.aru <- visitLimit
nsurveys.pc <- n_distinct(dfc$visit)
score <- samples$rebnut

y.aru <- matrix(NA, nsites, nsurveys.aru)

for (row in 1:nrow(MLcounts)) {
  y.aru[MLcounts$pointIndex[row], MLcounts$visit[row]] <- MLcounts$count[row]
}

y.pc <- matrix(NA, nsites, nsurveys.pc)

for (row in 1:nrow(yRBNU)) {
  y.pc[yRBNU$pointIndex[row], yRBNU$visit[row]] <- yRBNU$abun[row]
}


# define variables and make list of data ----------------------------------


data <- list(
  y.pc = y.pc,
  y.aru = y.aru,
  siteid = siteID,
  occid = occID,
  nsites = nsites,
  nsamples = nsamples,
  nsurveys.pc = nsurveys.pc,
  nsurveys.aru = nsurveys.aru,
  score = score
)


# Specify Model C in BUGS language
# Occupancy with false-positives, bioacoustics data with built-in
# species classification using Gaussian mixtures
cat(file = "modelC.txt", "
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
      y[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    for(j in 1:nsurveys.aru) { # Loop over occasions
      yARU[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
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

# Initial values
zst <- rep(1, nrow(y.pc))
gst <- sample(1:2, length(score), replace = TRUE)
gst[score > 0] <- 1
gst[score <= 0] <- 2
inits <- function() {
  list(
    mu = c(1, -1), sigma = 0.2, z = zst,
    psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
    lam = runif(1, 1, 2), ome = runif(1, 0, 0.4), g = gst
  )
}

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma", "Npos")

# MCMC settings
na <- 1000
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS (ART 1 min), gauge convergence and summarize posteriors
out3 <- jags(data, inits, params, "modelC.txt",
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
)


# visualize model results -------------------------------------------------

output <- as.data.frame(out3$summary[1:10, ])
