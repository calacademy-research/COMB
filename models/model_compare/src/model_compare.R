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
library(jagsUI)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_top1.R"))


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
threshold <- -3
aruVisitLimit <- 24

# Reading data
  
plan(multisession, workers = 20)

data <- future_map(speciesCodes, 
              ~ readCombined(
                  species = c(.x),
                  years = c(2021),
                  beginTime = dhours(6),
                  endTime = dhours(10),
                  visitLimit = aruVisitLimit,
                  visitAggregation = "file",
                  thresholdOptions = list(
                    value = threshold,
                    is.quantile = F
                  ),
                  squeeze = T
))

names(data) <- speciesCodes

# Slicing data
data_slice <- function(data_list, pc_days, aru_visits){
pc_days <- pc_days # number of pc days to subset
aru_visits <- aru_visits # number of ARU days to subset (possibly add amount of time per day to subset later)

pc_sub <- data_list[["indices"]][["visit.pc"]] %>% 
  filter(Visit <= pc_days)

aru_sub <- data_list[["indices"]][["visit.aru"]] %>% 
  filter(Visit_Index <= aru_visits)

data_list[["indices"]][["visit.aru"]] <- aru_sub

data_list[["indices"]][["visit.pc"]] <- pc_sub

# Adjusting the matrix indices
data_list[["y.ind"]] <- data_list[["y.ind"]][,1:pc_days, drop = F]
data_list[["y.aru"]] <- data_list[["y.aru"]][,1:aru_visits, drop = F]

# Adjusting the number of visits
data_list[["nsurveys.aru"]] <- aru_visits
data_list[["nsurveys.pc"]] <- pc_days


return(data_list)
}

# Defining all possible combinations of days and visits
# combinations <- expand_grid(speciesCodes, 1:3, 1:24) %>% 
#   rename("species" = 1, "pc" = 2, "aru" = 3)

combinations <- expand_grid(1:3, 1:24) %>% 
  rename("pc" = 1, "aru" = 2)

all_data <- pmap(list(combinations$pc, combinations$aru), 
                                ~ data_slice(data[["GCKI"]], ..1, ..2))

# all_data <- pmap(list(combinations$species, combinations$pc, combinations$aru), 
#               ~ data_slice(data[[..1]], ..2, ..3))

model <- function(data_list){
  data <- data_list
### BEGINNING OF MODEL COPY-PASTED
modelFile <- tempfile()
cat(file = modelFile, "
model {
  # Priors
 # psi ~ dunif(0, 1) # psi = Pr(Occupancy)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  r ~ dunif(1,100)  
  omega ~ dunif(0,15)
  lambda ~ dunif(0,25)
  beta0 ~ dnorm(0, 10) # Intercept for regression 
  beta1 ~ dnorm(0, 10) # Slope for regression 
  # Parameters of the observation model for the scores
  mu[1] ~ dunif(-1, 1.5)
  mu[2] ~ dunif(-3, 0)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])
  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*veg[i]
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    lam[i] <- lambda*z[i] + omega + 0.00001
    pnb[i] <- lam[i]/(r + lam[i])
    
    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    # GOF Point Count - Tukey-Freeman Discrepancy
    Tobs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    ySim[i] ~ dbin(p11 * z[i], n_v[i])
    Tsim0[i] <- (sqrt(ySim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*veg[i])/(1+exp(beta0 + beta1*veg[i]))
    zSim[i] ~ dbern(psiSim[i])
    
    # ARU - Negative binomial 
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dnegbin(pnb[i], r)  # n_surveys.aru = Total files processed
      
    #GOF ARU - Deviance
      LLobs0[i,j] <- logdensity.negbin(y.aru[i,j], pnb[i], r)
      y_aru_Sim[i,j] ~ dnegbin(pnb[i], r)
      LLsim0[i,j] <- logdensity.negbin(y_aru_Sim[i,j], pnb[i], r)
    }
    LLobs1[i] <- sum(LLobs0[i,])
    LLsim1[i] <- sum(LLsim0[i,])
    
  }
  # Likelihood part 2: feature score data
  for (i in 1:nsites) {
    site.prob[i] <- lambda*z[i]/(lambda*z[i]+omega) # Pr(sample is target species)
  }
  for(k in 1:nsamples) {
    # Sample specific covariate
    score[k] ~ dnorm(mu[g[k]], tau[g[k]]) # parameters are group specific
    probs[k,1] <- site.prob[siteid[k]]
    probs[k,2] <- 1 -site.prob[siteid[k]] # the prior class probabilities
    g[k] ~ dcat(probs[k,])
 #   N1[k] <- ifelse(g[k]==1, 1, 0)
  }
  # GOF assessment
  T_pc_obs <- sum(Tobs0)
  T_pc_sim <- sum(Tsim0)
  Dobs <- -2 * sum(LLobs1)
  Dsim <- -2 * sum(LLsim1)
  
  #GOF - Regression: Difference in mean vegetation
   cz1 <- sum(z)
   cz0 <- sum(1 - z)
   TvegObs <- sum(veg*z) / ifelse(cz1>0,cz1, 1)  - sum(veg*(1-z)) / ifelse(cz0>0,cz0, 1)
   czSim1 <- sum(zSim)
   czSim0 <- sum(1 - zSim)
   TvegSim <- sum(veg*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(veg*(1-zSim))/ifelse(czSim0>0,czSim0, 1)
  
  # Derived quantities
 # Npos <- sum(N1[])
}
")


# initialization
zst <- rep(1, data$nsites)
psit = runif(data$nsamples)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(1, -1), sigma = c(1, 1), z = zst,
    # psi = psit, 
    p11 = runif(1, 0.2, 0.8),
    lambda = runif(1, 1, 20), omega = runif(1, 0, 20), 
    g = gst, beta0 = 0, beta1 = 0
  )
}



# JAGS execution ----------------------------------------------------------

monitored <- c("beta0", "beta1", "p11", "lambda", "omega","r", "mu", "sigma", "T_pc_obs", "T_pc_sim", "Dobs", "Dsim",
               "TvegObs", "TvegSim")


# MCMC settings
na <- 1000
ni <- 8000
nt <- 1
nb <- 1000
nc <- 1

# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
n_v_per_site <- rowSums(!is.na(data$y.ind[, 1:data[["nsurveys.pc"]], drop = F])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:data[["nsurveys.pc"]], drop = F], na.rm = TRUE)
n_samp_per_site <- rowSums(!is.na(data$y.aru[, 1:data[["nsurveys.aru"]], drop = F])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:data[["nsurveys.aru"]], drop = F], na.rm = TRUE)


data2 <- append(data, list(n_v = n_v_per_site, y_pc_sum = y_pc_sum, n_s = n_samp_per_site, y_aru_sum = y_aru_sum))
# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))

covs <- read_csv(here("models/input/wide4havars.csv")) %>%
  mutate(Point = avian_point) %>% 
  dplyr::select(Point, mean_CanopyCover_2020_4ha) 

Veg <-left_join(as.data.frame(data$indices$point), covs, by = "Point" )
jagsData$veg <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))

set.seed(123)

jagsResult <- jags(jagsData, inits, monitored, modelFile,
                   n.adapt = na,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = F,
)

return(jagsResult)
### END OF MODEL COPY-PASTED
}

plan(multisession, workers = 10)

big_jags_out <- future_map(.x = 1:length(all_data), 
                    ~ model(all_data[[.x]]))

flatMeanAndRhat <- function(jagsResult) {
  means <- unlist(jagsResult$mean)
  rhats <- unlist(jagsResult$Rhat)
  names(rhats) <- as.character(
    map(names(rhats), 
        ~ paste("rhat", n, sep = "_")
    )
  )
  append(means, rhats) %>% 
    t() %>% 
    as.data.frame()
}

result_table <- map_df(1:length(all_data),
                       ~ flatMeanAndRhat(big_jags_out[[.x]])
)


