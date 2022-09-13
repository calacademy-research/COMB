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
source(here("models/src/model_read_lib_top3.R"))

# parameters --------------------------------------------------------------
# Species tried: HAWO, PIWO, BUSH, OSFL
speciesCode <- "BUSH" # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- -100
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
   # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    mu[i] <- lambda[i]*z[i] + omega + 0.00001
    z[i] ~ dbern(psi) # Latent occupancy states
    lambda[i] ~ dunif(0,25)
    
    # ARU
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dpois(mu[i])  # n_Total samples processed
    #GOF
      LLobs0[i,j] <- logdensity.pois(y.aru[i,j], (mu[i]))
      y_aru_Sim[i,j] ~ dpois(mu[i])
      LLsim0[i,j] <- logdensity.pois(y_aru_Sim[i,j], (mu[i]))
    }
    LLobs1[i] <- sum(LLobs0[i,])
    LLsim1[i] <- sum(LLsim0[i,])
  }
  # priors:
  psi ~ dunif(0, 1) # proportion of non-zeros
    
    omega ~ dunif(0,15)
  # GOF assessment
  Dobs <- -2 * sum(LLobs1)
  Dsim <- -2 * sum(LLsim1)
}")

# initialization
zst <- rep(1, data$nsites)
lambdat <- runif(data$nsites, 0, 25)
omegat <- runif(1, 0, 10)
inits <- function() {
  list(
    z = zst,
    psi = runif(1), 
    lambda = lambdat, 
    omegat = omegat
  )
}

n_samp_per_site <- rowSums(ifelse(!is.na(data$y.aru[, 1:24]),1,0)) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:24], na.rm = TRUE)
data2 <- append(data, list(n_v = n_v_per_site, y_pc_sum = y_pc_sum, n_s = n_samp_per_site, y_aru_sum = y_aru_sum))
jagsData <- within(data2, rm(indices))


monitored <- c("psi",  "lambda", "omega", "Dobs", "Dsim")

# MCMC settings
na <- 1000
ni <- 4000
nt <- 1
nb <- 1000
nc <- 6


jagsResult <- jags(jagsData, inits, monitored, modelFile,
                   n.adapt = na,
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
)

pp.check(jagsResult, "Dobs", "Dsim", main="Posterior Predictive Check\nDeviance - ARU")


##### Checking the ARU data itself - pay attention to: 
#  1) how many sites have 0/non-zero counts,
#  2) the mean (across files within sites)
#  3) the variance of counts (across files within sites) 
# - and if the variance is much larger than the mean

nonzero_mean <- function(x){
  return(mean(x[!(x==0) & !is.na(x)]))
}

nonzero_var <- function(x){
  return(var(x[!(x==0) & !is.na(x)]))
}

nz_count_stats <- cbind.data.frame(
  n_samp_per_site = n_samp_per_site,
  n_nz_count = rowSums(ifelse(!is.na(data$y.aru[, 1:24]) & 
                                     (!(data$y.aru[, 1:24])==0),1,0)),
  nz_mean = apply(data$y.aru[, 1:24], 1, nonzero_mean),
  nz_var = apply(data$y.aru[, 1:24], 1, nonzero_var)
)

print(nz_count_stats)

range = c(min(c(nz_count_stats$nz_mean, nz_count_stats$nz_var),na.rm = TRUE),
          max(c(nz_count_stats$nz_mean, nz_count_stats$nz_var),na.rm = TRUE))

ggplot(data = nz_count_stats, aes(x = nz_mean, y = nz_var)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name="Mean of non-zero counts", limits=range) +
  scale_y_continuous(name="Variance of non-zero counts", limits=range)

