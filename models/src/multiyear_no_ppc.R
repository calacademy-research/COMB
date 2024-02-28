# libraries ---------------------------------------------------------------
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)
library(abind)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated


# parameters --------------------------------------------------------------
speciesCode <- "BBWO" # must match prefiltering of dataML_model.csv
year <- c(2018,2019,2020,2021)
threshold <- 0
aruVisitLimit <- 32 # only consider this many ARU visits per site (ordered)


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
  years = year,
  beginTime = dhours(5),
  endTime = dhours(9),
  visitLimit = aruVisitLimit,
  visitAggregation = "file",
  thresholdOptions = list(value = threshold,
                          is.quantile = F),
  squeeze = T,
  logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
)
data$y.pc[1,1,,] # performs a quick check to see if pointcount data were read in properly
data$y.aru[1,,] # aru detections
data$score




# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  #p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) # commenting out p11 bc we are providing parameters for p below 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
  alpha0 ~ dunif(-5, 5)
  alpha1 ~ dunif(-5, 5)
  alpha2 ~ dunif(-5,5)
  alpha3 ~ dunif(-5, 5)
  alpha4 ~ dunif(-5,5)

   # Occupancy intercept on prob. scale
  for (t in 1:nyears){
    beta0[t] ~ dunif(-5, 5)
 }
 # beta0[2] ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
  # beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    for (t in 1:nyears){
      logit(psi[t,i]) <- beta0[t] + beta1*burn[i]*ifelse(t > 2, 1, 0)
      z[t,i] ~ dbern(psi[t,i]) # Latent occupancy states

      # Point count detection process
   
      for(j in 1:nsurveys.pc) {
        y.ind[t,i,j] ~ dbern(p11[t,i,j]*z[t,i]) # Note that the detection equation doesn't depend on year
        logit(p11[t,i,j]) <- alpha0 + alpha1*Time[t,i,j] + alpha2*Time[t,i,j]*Time[t,i,j] + alpha3*Date[t,i,j] + alpha4*Date[t,i,j]*Date[t,i,j]
      }
    
       # ARU - binomial
      p_aru[t,i] <- z[t,i]*p_aru11 + p_aru01 # Detection probability with false positives, no dependence on year
      for(j in 1:nsurveys.aru) {
        y.aru[t,i,j] ~ dbern(p_aru[t,i])  # n_surveys.aru = Total files processed
      }
    }
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[yearid[k],siteid[k]] + 1], tau[z[yearid[k],siteid[k]] + 1])
  }
    
    # derived quantity
 for (t in 1:nyears) {
  mean_psi <- mean(psi[t,])
}
")


# inits for FULL MODEL ----------------------------------------------------


zst <- matrix(1, data$nyears, data$nsites)
psit = matrix(runif(data$nyears*data$nsamples), data$nyears, data$nsites)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(-2, 0),
    sigma = c(1, 1),
    z = zst,
    # psi = psit,
    #p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    g = gst,
    beta0 = rep(0,data$nyears),
    beta1 = 0,
    reg_parm = 1
  )
}

monitored <- c(
  "mean_psi",
  "beta0",
  "beta1",
  #  "beta2",
  "alpha0",
  "alpha1",
  "alpha2",
  "alpha3",
  "alpha4",
  # "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma"
  
)

#monitored <- c("z")

# MCMC settings
na <- 100
ni <- 5000
nt <- 1
nb <- 1000
nc <- 6


jagsData <- within(data, rm(indices))

# SITE covars
site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point)

colnames(site_covars) # list of options for site-level covariates-- I chose measures of canopy cover and burn severity at the 4ha scale

covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

Veg <-left_join(as.data.frame(data$indices$point), covs, by = "Point") 

jagsData$burn <- Veg$mean_RAVGcbi_20182019_4ha

#VISIT covars


PC <- read_csv("./models/input/PC_delinted.csv")
PC$DateTime <- with_tz(PC$DateTime, tzone = "America/Los_Angeles")
PC$JDay <- as.numeric(format(PC$DateTime, "%j"))
PC$Time <- times(format(PC$DateTime, "%H:%M:%S"))
PC$JDay.s <- scale(PC$JDay)
PC$Time.s <- scale(PC$Time)

visit_pc <- PC %>% filter(year %in% c(2018,2019,2020,2021), birdCode_fk==speciesCode) %>% dplyr::select(year, point_ID_fk, visit, JDay, Time, JDay.s, Time.s)
pcTime <- as.matrix(visit_pc %>% group_by(year, point_ID_fk, visit) %>% dplyr::select(year, point_ID_fk, visit, Time.s) %>% pivot_wider(names_from = visit, values_from = Time.s))

pcDate <- as.matrix(visit_pc %>% group_by(year, point_ID_fk, visit) %>% dplyr::select(year, point_ID_fk, visit, JDay.s) %>% pivot_wider(names_from = visit, values_from = JDay.s))



jagsData$Time <-aperm(abind(pcTime[year==2018,3:5],pcTime[year==2019,3:5],pcTime[year==2020,3:5],pcTime[year==2021,3:5], along = 3), c(3,1,2))
jagsData$Time[is.na(jagsData$Time)] <- 0

jagsData$Date <- aperm(abind(pcDate[year==2018,3:5],pcDate[year==2019,3:5],pcDate[year==2020,3:5],pcDate[year==2021,3:5] ,along = 3), c(3,1,2))
jagsData$Date[is.na(jagsData$Date)] <- 0


set.seed(143617)

jagsResult_BBWO_multiyear <- jags(
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
