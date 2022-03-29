###########################################
##  These are some of the exercises from the Book,
##  Applied Heirarchical Modeling in Ecology, Vol 2
##  By Kery and Royle
##  Chapter 7. Exercises are cut and pasted here and most were debugged
##  To work on our systems
##  December, 2021
##  Requires the packages below and a separate installation of JAGS
###########################################



library(unmarked)
library(AHMbook)
library(ggplot2)
library(rjags)
library(jagsUI)


# Simulation settings set.seed(1)
# Initialize RNGs
nsites <- 200 # number of sites (i = 1, ..., nsites=M)
nsurveys <- 7 # number of visits (j = 1, ..., nsurveys=J)
psi <- 0.6 # expected occupancy probability
p <- 0.7 # detection probability (p_11)
fp <- 0.05 # false-positive error probability (p_10)

# Simulate occupancy states and encounter histories
z <- rbinom(nsites, 1, psi) # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for detections
for (i in 1:nsites) {
  pr_yequals1 <- p * z[i] + fp * (1 - z[i]) # p11 + p10
  y[i, ] <- rbinom(nsurveys, 1, pr_yequals1) # realized observations
}

# Number of false-positive detections per occasion
apply(y[z == 0, ] > 0, 2, sum)
# Number of false-negative detections per occasion
apply(y[z == 1, ] == 0, 2, sum)


type <- c(0, 7, 0)
# Build the unmarkedFrame
library(unmarked)
summary(umf <- unmarkedFrameOccuFP(y = y, type = type)) # not shown

largerp11 <- qlogis(c(0.5, 0.7, 0.1))
largerp10 <- qlogis(c(0.5, 0.1, 0.7))

(m1 <- occuFP(
  detformula = ~1, # model for p_11
  FPformula = ~1, # model for p_10
  stateformula = ~1, # model for psi
  data = umf, # umarkedFrameOccuFP object
  starts = largerp11
)) # add p_10 < p_11 constraint

(m2 <- occuFP(
  detformula = ~1, # model for p_11
  FPformula = ~1, # model for p_10
  stateformula = ~1, # model for psi
  data = umf, # umarkedFrameOccuFP object
  starts = largerp10
)) # add p_11 < p_10 constraint

cbind("m1" = plogis(coef(m1)), "m2" = plogis(coef(m2)))

# Look at the water vole data (in the AHMbook package)
data(waterVoles)
wv <- waterVoles
head(wv, 10)

# Make the false-positive umf (note not all sites surveyed in each year)
summary(wv.umf <- unmarkedFrameOccuFP(
  y = wv[, c("y1", "y2", "y3")],
  siteCovs = wv[, c("Year"), drop = FALSE], type = c(0, 3, 0)
)) # not shown

# 'Means' parameterization of these models
stvals <- list(
  "null" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.1)),
  "p11.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.1)),
  "p10.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.1, 0.1, 0.1)),
  "both.t" = qlogis(c(0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1))
)

cand.mods <- list(
  "null" = occuFP(
    detformula = ~1, FPformula = ~1,
    stateformula = ~ Year - 1, data = wv.umf, starts = stvals$null
  ),
  "p11.t" = occuFP(
    detformula = ~ Year - 1, FPformula = ~1,
    stateformula = ~ Year - 1, data = wv.umf, starts = stvals$p11.t
  ),
  "p10.t" = occuFP(
    detformula = ~1, FPformula = ~ Year - 1,
    stateformula = ~ Year - 1, data = wv.umf, starts = stvals$p10.t
  ),
  "both.t" = occuFP(
    detformula = ~ Year - 1, FPformula = ~ Year - 1,
    stateformula = ~ Year - 1, data = wv.umf, starts = stvals$both.t
  )
)

(modTab <- modSelFP(cand.mods))

# install.packages("AICcmodavg")
library(AICcmodavg)
aictab(cand.mods)


# Select the AIC-top model
topmod <- cand.mods$both.t

# Values for prediction
pred.df <- data.frame(Year = factor(c(2009, 2010, 2011)))

# Predict
p1 <- predict(topmod, type = "det", newdata = pred.df)
p2 <- predict(topmod, type = "fp", newdata = pred.df)
p3 <- predict(topmod, type = "state", newdata = pred.df)

# Create a data frame of predictions (code for plotting on book website)
preds <- rbind(p1, p2, p3)
preds$Type <- rep(factor(c("Detection", "False Positive", "Occupancy")), each = 3)
preds$Year <- rep(c("2009", "2010", "2011"), times = 3) # produce estimates for Fig. 7.2

ggplot(data = preds) +
  geom_point(mapping = aes(x = Year, y = Predicted, group = Type)) +
  facet_grid(. ~ Type)


# 7.3 mixtures of type I annd Type II data:

set.seed(129) # RNG seed
nsites <- 200 # number of sites
nsurveys1 <- 3 # number of occasions with Type 1 data
nsurveys2 <- 4 # number of occasions with Type 2 data
psi <- 0.6 # expected proportion of sites occupied
p <- c(0.7, 0.5) # detection prob of method 1 and method 2
fp <- 0.05 # false-positive error probability (p_10)

# Simulate latent occupancy states and data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
for (i in 1:nsites) {
  p1 <- p[1] * z[i] # certainly detection (method 1)
  p2 <- p[2] * z[i] + fp * (1 - z[i]) # uncertainly detection (method 2)
  y[i, 1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
  y[i, 4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data
}
# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 3), rep("2", 4)),
  nrow = nsites,
  ncol = nsurveys1 + nsurveys2, byrow = TRUE
)

type <- c(nsurveys1, nsurveys2, 0)

summary(umf1 <- unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method), type = type)) # not shown

(m1 <- occuFP(detformula = ~Method, FPformula = ~1, stateformula = ~1, data = umf1))


# Coefficients on the link (= "beta") scale
coef(m1)

# Coefficients on the probability (="real") scale
pred.df <- data.frame(Method = c("1", "2"))
round(rbind(
  "det" = predict(m1, type = "det", newdata = pred.df),
  "fp" = predict(m1, type = "fp", newdata = pred.df[1, , drop = F]),
  "state" = predict(m1, type = "state", newdata = pred.df[1, , drop = F])
), 3)

# 7.3 Joint modeling of type 1 and type 2 data()

# Load the data (from AHMbook) and summarize
data(EurasianLynx)
str(lynx <- EurasianLynx)

# Add the columns we need for analysis in unmarked
lynx$occ.1 <- 1
lynx$occ.2 <- 2
lynx$occ.3 <- 3
lynx$sYear <- standardize(lynx$Year)

# Extract the type 1 and type 2 data separately and bind them together
lynx1 <- lynx[lynx[, "type"] == "certain", ]
lynx2 <- lynx[lynx[, "type"] == "uncertain", ]
lynx <- cbind(lynx1[, c(2, 3:5)], lynx2[, 3:5])
colnames(lynx) <- c("site.nr", "y.1", "y.2", "y.3", "y.4", "y.5", "y.6")
occ <- cbind(
  lynx1[, c("occ.1", "occ.2", "occ.3")],
  lynx2[, c("occ.1", "occ.2", "occ.3")]
)
colnames(occ) <- c("occ.1", "occ.2", "occ.3", "occ.4", "occ.5", "occ.6")
lynx <- cbind(lynx, lynx1[, c("Year", "sYear", "Cntry")])


# Make the false-positive unmarkedFrame. Be sure to indicate type!
y <- lynx[, paste0("y.", 1:6)]
siteCovs <- lynx[, c("sYear", "Year", "Cntry")]
obsCovs <- list(occ = occ)
summary(lynx.umf <- unmarkedFrameOccuFP(
  y = y,
  siteCovs = siteCovs, obsCovs = obsCovs, type = c(3, 3, 0)
))

# Models to investigate trends in recording (takes approx. 15 mins)
cand.mods <- list(
  "p(c)fp(c)psi(c*t)" = occuFP(~Cntry, ~Cntry, ~1, ~ Cntry * sYear, data = lynx.umf),
  "p(c*t)fp(c)psi(c*t)" = occuFP(~ Cntry * sYear, ~Cntry, ~1, ~ Cntry * sYear, data = lynx.umf),
  "p(c)fp(c*t)psi(c*t)" = occuFP(~Cntry, ~ Cntry * sYear, ~1, ~ Cntry * sYear, data = lynx.umf),
  "p(c*t)fp(c*t)psi(c*t)" = occuFP(~ Cntry * sYear, ~ Cntry * sYear, ~1, ~ Cntry * sYear, data = lynx.umf),
  "p(c+t)fp(c+t)psi(c*t)" = occuFP(~ Cntry + sYear, ~ Cntry + sYear, ~1, ~ Cntry * sYear, data = lynx.umf),
  "p(c+t)fp(c+t)psi(c+t)" = occuFP(~ Cntry + sYear, ~ Cntry + sYear, ~1, ~ Cntry + sYear, data = lynx.umf),
  "p(c)fp(c)psi(c)" = occuFP(~Cntry, ~Cntry, ~1, ~Cntry, data = lynx.umf),
  "p(c*t)fp(c)psi(c)" = occuFP(~ Cntry * sYear, ~Cntry, ~1, ~Cntry, data = lynx.umf),
  "p(c)fp(c*t)psi(c)" = occuFP(~Cntry, ~ Cntry * sYear, ~1, ~Cntry, data = lynx.umf),
  "p(c*t)fp(c*t)psi(c)" = occuFP(~ Cntry * sYear, ~ Cntry * sYear, ~1, ~Cntry, data = lynx.umf),
  "p(c+t)fp(c+t)psi(c)" = occuFP(~ Cntry + sYear, ~ Cntry + sYear, ~1, ~Cntry, data = lynx.umf),
  "p(c)fp(c)psi(t)" = occuFP(~Cntry, ~Cntry, ~1, ~sYear, data = lynx.umf),
  "p(c*t)fp(c)psi(t)" = occuFP(~ Cntry * sYear, ~Cntry, ~1, ~sYear, data = lynx.umf),
  "p(c)fp(c*t)psi(t)" = occuFP(~Cntry, ~ Cntry * sYear, ~1, ~sYear, data = lynx.umf),
  "p(c*t)fp(c*t)psi(t)" = occuFP(~ Cntry * sYear, ~ Cntry * sYear, ~1, ~sYear, data = lynx.umf),
  "p(c+t)fp(c+t)psi(t)" = occuFP(~ Cntry + sYear, ~ Cntry + sYear, ~1, ~sYear, data = lynx.umf)
)

#  We compare the models using either the modSelFP or the aictab function.
(modTab <- modSelFP(cand.mods))

# Select and summarize the AIC-top model
(topmod <- cand.mods$`p(c+t)fp(c+t)psi(c*t)`)

# 7.4.1 MODELING CLASSIFIED FALSE POSITIVES IN unmarked
# Set parameter values of the simulation
set.seed(129) # RNG seed
nsites <- 200 # number of sites
nsurveys1 <- 3 # number of occasions with Type 1 data
nsurveys2 <- 4 # number of occasions with Type 2 data
psi <- 0.6 # expected proportion of are occupied
p <- c(0.7, 0.5) # detection prob of method 1 and method 2
fp <- 0.05 # false-positive error probability (p_10)
b <- 0.2 # probability y is recorded as certain

# Simulate the occupancy states and data
z <- rbinom(nsites, 1, psi)
y <- matrix(NA, nrow = nsites, ncol = nsurveys1 + nsurveys2)
for (i in 1:nsites) {
  p1 <- p[1] * z[i] # certainly detection (method 1)
  p2 <- p[2] * z[i] + fp * (1 - z[i]) # uncertainly detection (method 2)
  y[i, 1:3] <- rbinom(nsurveys1, 1, p1) # simulate method 1 data
  y[i, 4:7] <- rbinom(nsurveys2, 1, p2) # simulate method 2 data

  # Now introduce certain observations:
  pr.certain <- z[i] * y[i, 4:7] * b
  y[i, 4:7] <- y[i, 4:7] + rbinom(4, 1, pr.certain)
}

# Make a covariate to distinguish between the two methods
Method <- matrix(c(rep("1", 3), rep("2", 4)),
  nrow = nsites,
  ncol = nsurveys1 + nsurveys2,
  byrow = TRUE
)

#  Occasions 1e3 contain Type 1 and occasions 4e7 Type 3 data.
type <- c(nsurveys1, 0, nsurveys2)


# We create an unmarkedFrameOccuFP data object and fit the model. We have no covariates here
# and therefore no siteCovs or obsCovs either. With part of the data classified as certain, we can estimate the additional parameter b, and hence, now have four model formulas, each corresponding to a basic parameter: (1) detection probability, detformula (p11); (2) false-positive probability, FPformula (p10); (3) certainty probability, Bformula (b); and, finally, (4) occupancy probability, stateformula (j).
summary(umf2 <- unmarkedFrameOccuFP(y = y, obsCovs = list(Method = Method), type = type)) # not printed

(m2 <- occuFP(detformula = ~ -1 + Method, FPformula = ~1, Bformula = ~1, stateformula = ~1, data = umf2))

# We extract the coefficients and convert the estimates back to the natural scale:

# Coefficients on the link (= "beta") scale
coef(m2)

# Coefficients on the probability (="real") scale
pred.df <- data.frame(Method = c("1", "2"))
round(rbind(
  "det" = predict(m2, type = "det", newdata = pred.df),
  "fp" = predict(m2, type = "fp", newdata = pred.df[1, , drop = F]),
  "b" = predict(m2, type = "b", newdata = pred.df[1, , drop = F]),
  "state" = predict(m2, type = "state", newdata = pred.df[1, , drop = F])
), 3)

#  7.4.2 A GENERAL MULTITYPE MODEL WITH COVARIATES  - NOT added below...

# 7.5 BAYESIAN ANALYSIS OF MODELS WITH FALSE POSITIVES IN JAGS 421

# Random seed and simulation settings
set.seed(129, kind = "Mersenne")
nsites <- 200 # number of sites (i = 1, ..., M)
nsurveys <- 7 # number of visits (k = 1, ..., J)
psi <- 0.6 # expected psi
p <- 0.7 # detection probability (p_11)
fp <- 0.05 # false positive error probability (p_10)

# Simulate the latent states and the data
z <- matrix(NA, nrow = nsites, ncol = 1) # empty matrix for occ states
z[1:nsites] <- rbinom(nsites, 1, psi) # occupancy states
y <- matrix(NA, nrow = nsites, ncol = nsurveys) # empty matrix for det.
for (i in 1:nsites) {
  pr_yequals1 <- p * z[i] + fp * (1 - z[i]) # p11 + p10
  y[i, ] <- rbinom(nsurveys, 1, pr_yequals1) # realized observations
}

## Bundle data and summarize data bundle
str(bdata <- list(y = y, nsites = nrow(y), nsurveys = ncol(y)))

# Specify model in BUGS language
cat(file = "occufp.txt", "
model {

# Priors
psi ~ dunif(0, 1)
p ~ dunif(0, 1)
fp ~ dunif(0, 1)

# Likelihood and process model
for (i in 1:nsites) { z[i] ~ dbern(psi)
for (j in 1:nsurveys) {
  y[i,j] ~ dbern(z[i]*p + (1-z[i])*fp)
}
}
}
")

# We have to be careful to initialize the chains for the
# parameters in the correct region of the
# parameter space to avoid problems due to the multimodal likelihood.

# Initial values
zst <- apply(y, 1, max)
inits <- function() {
  list(z = zst, p = 0.7, fp = 0.05)
}

# Parameters monitored
params <- c("psi", "p", "fp")

# MCMC settings
na <- 1000
ni <- 5000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS (ART <1 min), assess convergence and summarize posteriors
library(jagsUI)
out1 <- jags(bdata, inits, params, "occufp.txt",
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
)
par(mfrow = c(2, 2))
traceplot(out1)
print(out1, 3)


# 7.6 MODELING FALSE POSITIVES FROM ACOUSTIC MONITORING DATA

# Simulation settings
set.seed(2019, kind = "Mersenne")
nsites <- 100 # Number of sites
nsurveys <- 5 # Number of replicates/occasions
psi <- 0.7 # Occupancy
p11 <- 0.5 # Detection probability at an occupied site
p10 <- 0.05 # False detection probability
lam <- 3 # Rate of true positives from ARU
ome <- 0.50 # Rate of false positives from ARU

# Simulate true occupancy states
z <- rbinom(nsites, 1, psi)

# Define detection probability
p <- z * p11 + (1 - z) * p10

# Simulate occupancy data and ARU count frequencies
yARU <- y <- K <- Q <- matrix(NA, nsites, nsurveys)
for (i in 1:nsites) {
  y[i, ] <- rbinom(nsurveys, 1, p[i]) # Detection/nondetection data
  K[i, ] <- rpois(nsurveys, lam * z[i]) # True positive detection frequency
  Q[i, ] <- rpois(nsurveys, ome) # False-positive detection frequency
  yARU[i, ] <- K[i, ] + Q[i, ] # Number of ARU detections
}

# Bundle and summarize data
str(bdata <- list(y = y, yARU = yARU, nsites = nsites, nsurveys = nsurveys))

# Specify Model A in BUGS language
cat(file = "modelA.txt", "
model {

# Priors
psi ~ dunif(0, 1)         # psi = Pr(Occupancy)
p10 ~ dunif(0, 1)         # p10 = Pr(y = 1 | z = 0)
p11 ~ dunif(0, 1)         # p11 = Pr(y = 1 | z = 1)
lam ~ dunif(0, 1000)
ome ~ dunif(0, 1000)

# Likelihood:process and observation models
for (i in 1:nsites) {
   z[i] ~ dbern(psi)                      # Occupancy status of site i
   p[i] <- z[i] * p11 + (1-z[i]) * p10    # false-positive detection model
   for(j in 1:nsurveys) {
      y[i,j] ~ dbern(p[i])                # Binary occupancy data
      yARU[i,j] ~ dpois(lam * z[i] + ome) # ARU detection frequency data
}
}
}
")

# Parameters monitored
params <- c("psi", "p10", "p11", "lam", "ome")

# MCMC settings
na <- 1000
ni <- 2000
nt <- 1
nb <- 1000
nc <- 3

# Call JAGS (tiny ART), gauge convergence and summarize posteriors
out1 <- jags(bdata, inits, params, "modelA.txt",
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
)
par(mfrow = c(2, 3))
traceplot(out1)
print(out1, 3)

# 7.6.2 MAKING USE OF ACOUSTIC VALIDATION DATA
# To integrate confirmation data, we first add some simulated
# confirmation data to our previously simulated data and then
# modify the model to accommodate the new validation data. For
# the simulation we use the valid_data function in AHMbook
# (modified from Chambert et al., 2017) which takes the observed
# ARU counts and requires the number of true positives to be input,
# from which some are (randomly) validated...

# The function is applied as follows to our simulated data set based on J (== nsurveys) replicates:
k <- n <- matrix(NA, nrow = nsites, ncol = nsurveys)
for (j in 1:nsurveys) {
  bleen <- valid_data(yARU[, j], tp = K[, j], n.valid = 50)
  n[, j] <- bleen$n
  k[, j] <- bleen$k
}
# The output objects are the number of samples checked for each site (n)
# and the number of confirmed positives (k)

# 7.6.3 THE COUPLED CLASSIFICATION MODEL: A NEW MODEL FOR INTEGRATING SPECIES CLASSIFICATION FROM BIOACOUSTICS DATA WITH OCCUPANCY MODELS

# Simulation settings
set.seed(2019)
nsites <- 100 # Number of sites
nsurveys <- 5 # Number of replicates
psi <- 0.7 # Occupancy
p11 <- 0.5 # Detection probability at an occupied site
p10 <- 0.05 # False-positive rate for occupancy surveys
lam <- 1 # Rate at which target species vocalizations are obtained
ome <- 3 # Rate at which non-target species are obtained
mu.mix.target <- 1 # mean of feature for target species
mu.mix.nontarget <- -1 # ... for non-target species
sd.mix <- 0.5 # sd of feature for both

# Simulate true occupancy states
z <- rbinom(nsites, 1, psi)

# Define detection probability, with false positives
p <- z * p11 + (1 - z) * p10

# Simulate occupancy data and ARU count frequencies
# Main simulation loop
covdata <- NULL
yARU <- y <- K <- Q <- matrix(NA, nsites, nsurveys)
for (i in 1:nsites) {
  for (j in 1:nsurveys) {
    y[i, j] <- rbinom(1, 1, p[i]) # Detection/nondetection data
    K[i, j] <- rpois(1, lam * z[i]) # True detection frequency
    Q[i, j] <- rpois(1, ome) # False-positive detection frequency
    yARU[i, j] <- K[i, j] + Q[i, j] # Number of ARU detections

    if (yARU[i, j] > 0) { # now simulate some scores and store them for output
      K.score <- rnorm(K[i, j], mu.mix.target, sd.mix)
      Q.score <- rnorm(Q[i, j], mu.mix.nontarget, sd.mix)
      covdata <- rbind(covdata, cbind(
        rep(i, yARU[i, j]), rep(j, yARU[i, j]),
        c(K.score, Q.score), c(rep(1, K[i, j]), rep(2, Q[i, j]))
      ))
    }
  }
}
colnames(covdata) <- c("siteID", "occID", "score", "truth")

# Harvest the data
siteid <- covdata[, 1]
occid <- covdata[, 2]
score <- covdata[, 3]
truth <- covdata[, 4]

# Bundle and summarize data
str(bdata <- list(
  y = y, score = score, yARU = yARU, siteid = siteid,
  occid = occid, nsamples = length(score), nsites = nsites,
  nsurveys = nsurveys
))

# Specify Model C in BUGS language
# Occupancy with false-positives, bioacoustics data with built-in
# species classification using Gaussian mixtures
cat(file = "modelC.txt", "
model {
# Priors
psi ~ dunif(0, 1)                # psi = Pr(Occupancy)
p10 ~ dunif(0, 1)                # p10 = Pr(y = 1 | z = 0)
p11 ~ dunif(0, 1)                # p11 = Pr(y = 1 | z = 1)
lam ~ dunif(0, 1000)             # lambda: rate of target-species calls detected
ome ~ dunif(0, 1000)             # omega: rate of non-target detections

# Parameters of the observation model for the scores
mu[1] ~ dnorm(0, 0.01)
mu[2] ~ dnorm(0, 0.01)
sigma ~ dunif(0, 10)
tau <- 1 / (sigma * sigma)

# Likelihood part 1: detection data and ARU counts
for (i in 1:nsites) {                     # Loop over sites
  z[i] ~ dbern(psi)                       # Latent occupancy states
  p[i] <- z[i]*p11 + (1-z[i])*p10         # Detection probability
  site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species)
  for(j in 1:nsurveys) {                  # Loop over occasions
    y[i,j] ~ dbern(p[i])                  # Observed occ. data (if available)
    yARU[i,j] ~ dpois(lam*z[i] + ome)     # Total samples processed
    }
}

# Likelihood part 2: feature score data
for(k in 1:nsamples) {                       # Sample specific covariate
  score[k] ~ dnorm(mu[g[k]], tau)            # parameters are group specific
  probs[k,1] <- site.prob[siteid[k]]
  probs[k,2] <- 1 - site.prob[siteid[k]]     # the prior class probabilities
  g[k] ~ dcat(probs[k,])
  N1[k] <- ifelse(g[k]==1, 1, 0)
}

# Derived quantities
Npos <- sum(N1[])
}
")

# Initial values
zst <- rep(1, nrow(y))
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
out3 <- jags(bdata, inits, params, "modelC.txt",
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE
)
par(mfrow = c(2, 3))
traceplot(out3)
print(out3, 3)
