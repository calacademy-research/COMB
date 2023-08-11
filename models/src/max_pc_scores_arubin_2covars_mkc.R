# includes simulations for site covariates
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
source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated

# parameters --------------------------------------------------------------
speciesCode <- "PIWO" # must match prefiltering of dataML_model.csv
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

# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  beta0 ~ dunif(-5, 5) # Intercept for occupancy logistic regression
  beta1 ~ dunif(-5, 5)
  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*cover[i] + beta2*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states

    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    # GOF Point Count - Tukey-Freeman Discrepancy
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*cover[i] + beta2*burn[i])/(1+exp(beta0 + beta1*cover[i] + beta2*burn[i])) 
    zSim[i] ~ dbern(psiSim[i])

    # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
    }
    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  # GOF assessment
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
  TvegObs <- sum(cover*z) / ifelse(cz1>0,cz1, 1)  - sum(cover*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(cover*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(cover*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites
  mean_psi <- mean(psi)
 
  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta2 * Xburn[k] # psi predictions for severity, holding canopy cover at mean
    logit(psi.pred.cov.ub[k]) <- beta0 + beta1 * Xcover[k] # psi predictions across cover, holding burn = 0 (unburned)
    logit(psi.pred.cov.b[k]) <- beta0 + beta1 * Xcover[k] + beta2 * 1 # psi predictions across cover, holding burn = 1 (avg burn)
    
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

# TODO: tighten up this code; someday automatically do this with model_read_lib or similar script? 
# TODO: ...and also include visit-level covariates?

site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point)

colnames(site_covars) # list of options for site-level covariates-- I chose measures of canopy cover and burn severity at the 4ha scale

site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point) # filter to only points with avian data

colnames(site_covars) # list of options for site-level covariates to choose from-- with select() below I chose measures of canopy cover and burn severity at the 4ha scale
covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

Veg <-left_join(as.data.frame(data$indices$point), covs, by = "Point") # ensures veg vars are in the same point order as the count data
Veg$cover.s <- scale(Veg$mean_CanopyCover_2020_4ha)
Veg$lgTree.s <- scale(Veg$Perc_LTg22mHt_2020_4ha)

jagsData$cover <- as.numeric(Veg$cover.s) # needs as.numeric() bc the vector in veg retains scaling information 
jagsData$burn <- as.numeric(Veg$mean_RAVGcbi_20182019_4ha) 
jagsData$lgTree <- as.numeric(Veg$lgTree.s)

summary(jagsData$burn)
summary(jagsData$cover)
jagsData$Xburn = seq(min(jagsData$burn), max(jagsData$burn), length.out=100)
jagsData$Xcover = seq(min(jagsData$cover), max(jagsData$cover), length.out=100)

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
  "beta1",
  "beta2",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu", # score distributions 
  "sigma", # sd for score distributions
  "NOcc", "PropOcc",
  "psi.pred.burn",
  "psi.pred.cov.ub",
  "psi.pred.cov.b",
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

## Quick QC checks
print(jagsResult)
# table(apply(data$y.pc[, 1, , ], 1, max), apply(data$y.aru, 1, max)) # number of sites with detections by ARU versus point count
# (site_pos <-
#     cbind(apply(data$y.pc[, 1, , ], 1, max), apply(data$y.aru, 1, max)))
# mean(apply(site_pos, 1, max), na.rm = TRUE) # average number of sites with positive detections from either source
# plogis(jagsResult$mean$beta0) # Mean of the posterior of the mean occupancy rate psi across sites from model
# jagsResult$mean$mean.psi
# 
# mean(apply(data$y.pc[, 1, , ], 1, max), na.rm = TRUE)
# mean(apply(data$y.aru, 1, max), na.rm = TRUE)

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


b1 <- as.data.frame(jagsResult$summary[13:112,]) # predicted psi for values of Xburn
b1$Xburn <- jagsData$Xburn
ggplot(b1) +
  geom_ribbon(aes(x=Xburn, ymin=`2.5%`, ymax=`97.5%`, alpha=0.2)) +
  geom_line(aes(x=Xburn, y=mean)) +
  #facet_wrap(~spp) +
  theme_classic() +
  labs(y="predicted occupancy probability (psi)", x = "RAVG score", title="Black-backed Woodpecker") +
  theme(legend.position = "none")

jagsData$Xcover

c_unb <- as.data.frame(jagsResult$summary[113:212,]) # predicted psi for values of Xburn
c_unb$Xburn <- jagsData$Xcover
ggplot(c_unb) +
  geom_ribbon(aes(x=Xburn, ymin=`2.5%`, ymax=`97.5%`, alpha=0.2)) +
  geom_line(aes(x=Xburn, y=mean)) +
  #facet_wrap(~spp) +
  theme_classic() +
  labs(y="predicted occupancy probability (psi)", x = "Cover", title="Black-backed Woodpecker") +
  theme(legend.position = "none")

c_b <- as.data.frame(jagsResult$summary[213:312,]) # predicted psi for values of Xburn
c_b$Xburn <- jagsData$Xcover
ggplot(c_b) +
  geom_ribbon(aes(x=Xburn, ymin=`2.5%`, ymax=`97.5%`, alpha=0.2)) +
  geom_line(aes(x=Xburn, y=mean)) +
  #facet_wrap(~spp) +
  theme_classic() +
  labs(y="predicted occupancy probability (psi)", x = "Cover", title="Black-backed Woodpecker") +
  theme(legend.position = "none")
