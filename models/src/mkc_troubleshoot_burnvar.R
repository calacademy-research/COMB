# binary fire variable with some dummy tests

# libraries ---------------------------------------------------------------
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated

# parameters --------------------------------------------------------------
speciesCode <- "BBWO" # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- 0
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered)


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
  years = c(year),
  beginTime = dhours(6),
  endTime = dhours(10),
  visitLimit = aruVisitLimit,
  visitAggregation = "file",
  thresholdOptions = list(value = threshold,
                          is.quantile = F),
  squeeze = T,
  logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
)
data$y.pc[,1,,] # quick look at point count data
data$y.aru # aru detection data
data$score # score data

site_pos <- data.frame(cbind(apply(data$y.pc[, 1, , ], 1, max, na.rm = T), apply(data$y.aru, 1, max)))
colSums(site_pos, na.rm=T) # 12 points in PC data, 27 points in ARU data
colnames(site_pos) <- c("naivePC", "naiveARU")
site_pos$Point <- as.numeric(rownames(site_pos))

# simple GLMs -------------------------------------------------------------

# first let's try a simple glm to model a single covariate-- burn severity-- on the pc data, all visits collapsed

site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point) # filter to only points with avian data

colnames(site_covars) # list of options for site-level covariates to choose from-- with select() below I chose measures of canopy cover and burn severity at the 4ha scale
covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point) %>% 
  #mutate(CanopyCover = mean_CanopyCover_2020_4ha/100) %>% 
  rename(LgTree = Perc_LTg22mHt_2020_4ha, CBI = mean_RAVGcbi_20182019_4ha, CanopyCover = mean_CanopyCover_2020_4ha) %>% mutate(CanopyCover = CanopyCover/100)

pcvis <- data$y.pc[,1,,]

colnames(pcvis) <- c("v1", "v2", "v3")
d <- cbind(pcvis, covs) %>% select(Point, v1, v2, v3, CBI) %>%
  pivot_longer(cols = v1:v3)
d$value[is.na(d$value)] <- 0 # removing NA's now just for simplicity

d <- left_join(site_pos, covs)

summary(glm1 <- glm(naivePC ~ CBI, data = d, family="binomial")) # definitely strong relationship
sim <- data.frame(CBI = seq(0, max(d$CBI), len=200))
sim$dpred <- predict(glm1, sim, type = "response")
plot(naivePC~CBI, d)
lines(dpred ~ CBI, sim)
# now try with ARU data

summary(glm2 <- glm(naiveARU ~ CBI, data = d, family="binomial")) # no relationship

sim$dpred2 <- predict(glm2, sim, type = "response")
plot(naiveARU ~ CBI, d, col="red4", pch=2) # aru raw data
points(naivePC+0.02 ~ CBI, d, pch=1, col="blue4") # pc raw data
lines(dpred2 ~ CBI, sim, col="red4") # prediction for aru data
lines(dpred ~ CBI, sim, col="blue4") # prediction for pc data 
legend(x = 1.7, y=0.9,
       legend = c("naive PC", "naive ARU"),
       pch = c(1, 2),
       col = c("blue4", "red4"))

# so, just looking at raw data (not accounting for false + or false -!), the ARU and PC data predict different presence probabilities. The slope is much steeper on the PC data, which is what I'd expect given that BBWO are rarely detected by human observers outside of burn areas, and ARUs may have some false +s. But there could be underlying biology here, not just false positives. Are ARUs picking up BBWO calls in green forest? If false + error is uniform across the landscape (which I assume it should be), shouldn't the slopes be more or less equal?

# next, let's try with the (continuous) burn variable as the only covariate and see what happens.

# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) for PCs
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) for ARUs
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p01 = Pr(y = 1 | z = 0)
  
# NEW! I changed beta0 and beta1 priors from dnorm(0,10) to dnorm(0,0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 0.01) # changed to increase variance
  mu[2] ~ dnorm(-2, 0.01)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states

    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
    }
 }
  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }
  
# I moved all the GOF assessment to its own section so I could read both it and the main model more easily
  # GOF assessment
  
  for (i in 1:nsites) {
      # GOF Point Count - Tukey-Freeman Discrepancy
        T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
        y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
        T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

      # GOF - Regression: Simulations
        psiSim[i] <- exp(beta0 + beta1*burn[i])/(1+exp(beta0 + beta1*burn[i])) 
        zSim[i] ~ dbern(psiSim[i])
        
      # GOF ARU Count - Tukey-Freeman Discrepancy
        T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
        y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
        T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }
        
      # GOF for ARU scores
  for (k in 1:nsamples) {
        LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
        score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
        LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }
  
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  mean_psi <- mean(psi)
 
  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for severity
  }
}
")


# Add GoF inputs to data list -----------------------------------------------

n_v_per_site <-
  rowSums(!is.na(data$y.ind[, 1:3])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
n_samp_per_site <-
  rowSums(!is.na(data$y.aru[, 1:24])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:24], na.rm = TRUE) # number of samples where y_aru=1

data2 <-
  append(
    data,
    list(
      n_v = n_v_per_site,
      y_pc_sum = y_pc_sum,
      n_s = n_samp_per_site,
      y_aru_sum = y_aru_sum
    )
  )


# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))


# Add site-level covariates to data list ----------------------------------

covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

Veg <-left_join(as.data.frame(data2$indices$point), covs, by = "Point") # ensures veg vars are in the same point order as the count data
#Veg$cover.s <- scale(Veg$mean_CanopyCover_2020_4ha)
#Veg$lgTree.s <- scale(Veg$Perc_LTg22mHt_2020_4ha)

#jagsData$cover <- as.numeric(Veg$cover.s) # needs as.numeric() bc the vector in veg retains scaling information 
jagsData$burn <- as.numeric(Veg$mean_RAVGcbi_20182019_4ha) 

summary(jagsData$burn)
#summary(jagsData$cover)
jagsData$Xburn = seq(0, 2.5, length.out=100)
#jagsData$Xcover = seq(-2.75, 1.8, length.out=100)

# Set monitored parameters and starting values ----------------------------
# initialization
zst <- rep(1, data$nsites)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(-2, 0),
    sigma = c(1, 1),
    z = zst,
    mean.psi = runif(1),
    p11 = runif(1, 0.2, 0.8),
    g = gst,
    reg_parm = 1
  )
}
monitored <- c(
  "beta0",
  "beta1", # continuous burn variable
#  "beta2",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu", # score distributions 
  "sigma", # sd for score distributions
  "NOcc", "PropOcc",
  "psi.pred.burn",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "TvegObs",
  "TvegSim",
  "Dobs",
  "Dsim"
)

# MCMC settings
na <- 1000
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6

# JAGS execution ----------------------------------------------------------

#set.seed(123)

jagsResult <- jags(
  jagsData,
  inits,
  monitored,
  modelFile,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)

print(jagsResult)
plot(jagsResult)


# Posterior Predictive Checks ---------------------------------------------

pp.check(
  jagsResult,
  "TvegObs",
  "TvegSim",
  main =
    paste(
      "Posterior Predictive Check\n
       Logistic Regression: Difference in Occupied Means\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult,
  "T_pc_obs",
  "T_pc_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult,
  "T_aru_obs",
  "T_aru_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult,
  "Dobs",
  "Dsim",
  main =
    paste(
      "Posterior Predictive Check\nDeviance - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

b1 <- as.data.frame(jagsResult$summary[12:111,]) # predicted psi for values of Xburn
b1$Xburn <- jagsData$Xburn
ggplot(b1) +
  geom_ribbon(aes(x=Xburn, ymin=`2.5%`, ymax=`97.5%`, alpha=0.2)) +
  geom_line(aes(x=Xburn, y=mean)) +
  #facet_wrap(~spp) +
  theme_classic() +
  labs(y="predicted occupancy probability (psi)", x = "RAVG score", title="Black-backed Woodpecker") +
  theme(legend.position = "none")


### so-- with the old priors on beta0 and beta1 (dnorm(0,10), the model converged well, but did not predict a strong association between burn severity and occupancy. BUT, with dunif(0.01), it also converged well and does find a strong relationship. Also, the first posterior predictive check (the one that checks the burn variable) with the first priors were failing/close to it, while the PPCs with the second priors are much better. The PPC on the ARU data is exactly 1 or very close to it... don't we want that value to be close to neither 0 nor 1?

# I think the wide CI at the higher end of the scale probably has to do with how relatively few severely burned sites there are compared to sites with RAVG score = 0 or very low. Could try with burn as a categorial variable


# burn as categorical -----------------------------------------------------
###### ...with code as-is
###  first I'm going to try just changing the input itself without changing the code
jagsData$burn <- ifelse(jagsData$burn > 0, 1, 0)
# go back and run the jags code

jagsResult2 <- jags(
  jagsData,
  inits,
  monitored,
  modelFile,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)

print(jagsResult2)

# this model failed its first PPC (difference in occupied means)
b2 <- as.data.frame(jagsResult2$summary[12:111,])



# binary burn diff intercepts ---------------------------------------------

# have to rewrite JAGS code - this currently overwrites the previous model. also I didn't recode the PPCs for this one yet.

modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) for PCs
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) for ARUs
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p01 = Pr(y = 1 | z = 0)
 
 # NEW: new priors for (binary) burn variable
  for(b in 1:2) {
  pOcc[b] ~ dbeta(1,1)
  bBurn[b] <- logit(pOcc[b])
  }

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- bBurn[burn[i]]
    z[i] ~ dbern(psi[i]) # Latent occupancy states

    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
    }
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  mean_psi <- mean(psi)
  burnfx <- bBurn[2] - bBurn[1]

}
")

jagsData$burn <- jagsData$burn+1 # make it 1 and 2 to match the code

monitored <- c(
  "bBurn",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu", # score distributions 
  "sigma", # sd for score distributions
  "NOcc", "PropOcc",
  "burnfx"
)

# MCMC settings
na <- 1000
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6

# JAGS execution ----------------------------------------------------------

#set.seed(123)

jagsResult3 <- jags(
  jagsData,
  inits,
  monitored,
  modelFile,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)

## Quick QC checks
print(jagsResult3)

bin_out <- as.data.frame(jagsResult3$summary)

plogis(-2.123)
plogis(-0.079)

bin_sims <- as.data.frame(jagsResult3$sims.list[1:2])

boxplot(plogis(bin_sims$bBurn.1), plogis(bin_sims$bBurn.2))

        