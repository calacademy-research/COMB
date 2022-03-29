# Avian Community Occupancy Model for Caples Creek watershed BEFORE prescribed fire
# This code prepares the data for a model without covariates in the occupancy structure for comparison with models containing veg structure covariates. The detection process is modeled by the linear and quadratic effects of date

# NOT RUNNING as of 2022-03-17: matrix creation DONE (using reshape2); multiyear JAGS code NOT working

# in the future, will have separate scripts for single-year, multi-year, before-after models, all of which should run autonomously, generate separate jags files in models/jags, and uniquely named outputs in models/output

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(chron)
library(jagsUI)
library(coda)
library(lubridate)
library(googledrive)
library(here)

# Read in data ------------------------------------------------------------

if (dir.exists(here("models/jags/")) == FALSE) {
  dir.create(here("models/jags"))
}
if (dir.exists(here("models/output/")) == FALSE) {
  dir.create(here("models/output"))
}

# TODO []: link to Google Drive &/or a link to `read_PC.R`. The Google Drive code in readPC.R does not work when using non-CalAcademy login. For now, doing a hard download.
source(here("caples_functions.R"))
drive_auth()
#
drive_sync(here("point_counts/data_ingest/input/"), "https://drive.google.com/drive/u/0/folders/13bXkNCZzFwGC8H4k-CJ4Cf3tYJ3t18zW")
#
# # Reading inputs -------------------
PointC <- fread(here("point_counts/data_ingest/output/PointC_2022-03-09.csv"))

# source(here("point_counts/data_ingest/src/readPC.R")) # trying to get this code to run that script

# REMOVE double-observer counts and other errors/duplicates, filter data by distance, etc  ---------------

# TODO []: decide how to choose surveys if >3 (using weather data? etc)
# list of how many surveys per point per season in the dataset
PCTC %>%
  select(point_ID_fk, Date) %>%
  mutate(Year = lubridate::year(mdy(Date))) %>%
  filter(Year != 2017) %>%
  distinct() %>%
  group_by(Year, point_ID_fk) %>%
  tally() %>%
  spread(Year, n) %>%
  View()

# find surveys with duplicate entries and investigate cause (double-observer, error in FMPro entry, etc)

PC <- PointC %>%
  filter(
    dist_detect != "flyover", # comment these on/off to restrict distances
    dist_detect != "100m+",
    dist_detect != "50to100m",
    point_ID_fk != "1072", # missing veg data
    birdCode_fk != "UNKN", birdCode_fk != "0", !is.na(birdCode_fk),
    birdCode_fk != "XXHU", birdCode_fk != "XXWO", # comment on/off to include XX__ entries
    year(DateTime) != 2017,
    observer_fk != "IASH", observer_fk != "MASC"
  ) %>%
  mutate(year = year(DateTime)) %>%
  group_by(pointCount_ID_fk, point_ID_fk, DateTime, year, observer_fk, birdCode_fk) %>%
  summarise(abun = n()) %>%
  spread(birdCode_fk, value = abun, fill = 0) %>%
  pivot_longer(AMDI:WISA, names_to = "birdCode_fk", values_to = "abun")

visits <- PointC %>%
  mutate(year = year(DateTime)) %>%
  filter(
    point_ID_fk != "1072",
    observer_fk != "IASH", observer_fk != "MASC",
    birdCode_fk != "UNKN", birdCode_fk != "XXHU", birdCode_fk != "XXWO", !is.na(birdCode_fk),
    year != 2017
  ) %>%
  group_by(point_ID_fk, year, DateTime) %>%
  summarise(rich = n_distinct(fbirdName)) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(visit = seq_along(DateTime))

extraVisits <- visits %>% filter(visit > 3)

dfc <- full_join(PC, visits) %>%
  ungroup() %>%
  dplyr::select(abun, observer_fk, DateTime, point_ID_fk, year, birdCode_fk, visit) %>%
  filter(visit < 5)

dups <- dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit) %>%
  summarise(nspec = n_distinct(birdCode_fk))

dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit, DateTime) %>%
  arrange(DateTime) %>%
  summarise(nspec = n_distinct(birdCode_fk)) %>%
  arrange(DateTime) # all of these are double-observer/training counts

dfc <- dfc[which(duplicated(dfc[, 4:7]) == FALSE), ]

# check this is 0
dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit) %>%
  summarise(nspec = n_distinct(birdCode_fk)) # okay

dfc <- dfc %>% filter(!is.na(birdCode_fk))

dfc$observer_fk[dfc$observer_fk == "MC"] <- "MKC"
unique(dfc$observer_fk)

# visitswide <- dfc %>% dplyr::select(-observer_fk, -DateTime) %>% pivot_wider(names_from=visit, values_from = abun)

# visitswide %>% filter(!is.na(`4`)) %>% group_by(point_ID_fk, year) %>% summarise(n_distinct(birdCode_fk))



# (4) Format data in a spreadsheet format into a 4D array # Determine required dimensions of 4D array
nsite <- length(unique(dfc$point_ID_fk))
nvisit <- length(unique(dfc$visit))
nyear <- length(unique(dfc$year))
nspec <- length(unique(dfc$birdCode_fk))
# Prepare array and pre-fill array with NAs
# BROKEN HERE! TODO []: below code from AHM2 to create the 4D array needed won't work

# CaplesArray <- array(NA, dim = c(nsite, nvisit, nyear, nspec))

# error message for below for loop: "Error in `[<-`(`*tmp*`, tmp$point_ID_fk[i], tmp$visit[i], tmp$year[i],  : subscript out of bounds"
# for(i in 1:nrow(tmp)) {
#  CaplesArray[tmp$point_ID_fk[i], tmp$visit[i], tmp$year[i], tmp$birdCode_fk[i]] <- tmp$abun[i] }


# PCmat currently melts to y using reshape, but not in the base code above
# potentially important data: visits >3? weather data as detection variables?

library(reshape2)

# PCmat %>% dplyr::select(point_ID_fk, year, visit, birdCode_fk, abun)

y <- melt(dfc, id.var = c("birdCode_fk", "point_ID_fk", "visit", "year"), measure.var = "abun") %>% acast(point_ID_fk ~ visit ~ year ~ birdCode_fk)

# Covariates

# Visit covariates (for detection model)

# DATE matrix

dates1 <- visits %>%
  dplyr::select(-rich) %>%
  filter(visit < 5) %>%
  mutate(JDay = as.numeric(format(DateTime, "%j"))) %>%
  select(-DateTime) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "JDay") %>%
  acast(point_ID_fk ~ visit ~ year)

# standardize dates
mdate <- mean(dates1, na.rm = TRUE)
sddate <- sqrt(var(dates1[1:length(dates1)], na.rm = TRUE))
date1 <- (dates1 - mdate) / sddate
date2 <- date1 * date1
# date1 <- as.matrix(date1)
# date2 <- as.matrix(date2)

date1[is.na(date1)] <- 0
date2[is.na(date2)] <- 0

# TIME matrix [under construction - issue with time zone and converting to list of matrices]

times <- visits %>%
  dplyr::select(-rich) %>%
  filter(visit < 5) %>%
  mutate(Time = times(format(DateTime, "%H:%M:%S"))) %>%
  select(-DateTime) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "Time") %>%
  acast(point_ID_fk ~ visit ~ year)
# times1 <- as.matrix(times[,3:length(times)])

# standardize times
mtime <- mean(times, na.rm = TRUE)
sdtime <- sqrt(var(times[1:length(times)], na.rm = TRUE))
time1 <- (times - mtime) / sdtime
time2 <- time1 * time1
# time1 <- as.matrix(time1)
# time2 <- as.matrix(time2)

time1[is.na(time1)] <- 0
time2[is.na(time2)] <- 0 # for now, setting these to 0

# vector for number of visits per site j
nvisit <- visits %>%
  filter(visit < 5) %>%
  group_by(point_ID_fk, year) %>%
  summarise(nvisits = max(visit))
K <- as_vector(nvisit[, 3])
# K <- 3


# Site Covariates ---------------------------------------------------------

# the following two csvs were saved from the 3/15/2021 version of 'Caples Points Metadata.xlsx' from shared GoogleDrive-- in the future could add the googledrive package to read in directly?
env <- read_csv("data/caples_env_1ha_20210315.csv")
utms <- read_csv("data/Wildlife_Sampling_Points.csv")
head(env)

prefire <- env %>%
  dplyr::select(point_d, mean_Cover_2018_1ha, mean_Height_2018_1ha, mean_NBR_2018_1ha, Perc_NtoE_2020_1ha, Binary_NtoE_2020_1ha, mean_Elevation_2020_1ha) %>%
  left_join(utms, by = c("point_d" = "point_id")) %>%
  left_join(nvisit, by = c("point_d" = "point_ID_fk")) %>%
  arrange(point_d)

siteList <- unique(PointC$point_ID_fk)

setdiff(prefire$point_d, siteList) # these three aren't point count locations
setdiff(siteList, prefire$point_d) # 1072 is missing in the metadata for some reason

siteData <- prefire %>%
  filter(point_d %in% unique(visits$point_ID_fk)) %>%
  separate(SZ_DNS2, c("Size", "Density"), sep = 1)

J <- n_distinct(siteData$point_d)

siteData %>% ggplot() +
  geom_boxplot(aes(x = Size, y = mean_Height_2018_1ha)) +
  labs(main = "Height and Size", y = "mean tree height (1ha)", x = "Tree Size Class (Becky's data)")
siteData %>% ggplot() +
  geom_boxplot(aes(x = Density, y = mean_Cover_2018_1ha)) +
  labs(mean = "Density and Cover", y = "mean cover (1ha)", x = "Tree Density Class (Becky's data)")
siteData %>% ggplot() +
  geom_boxplot(aes(x = Density, y = mean_NBR_2018_1ha)) +
  labs(main = "Density and NBR", y = "NBR (1 ha)", x = "Tree Density Class (Becky's data)")

# Cover <- scale(siteData$mean_Cover_2018_1ha)
Height <- scale(siteData$mean_Height_2018_1ha)
# NBR <- scale(siteData$mean_NBR_2018_1ha)
El <- scale(siteData$mean_Elevation_2020_1ha)


# Build matrix for JAGS ---------------------------------------------------



# occupancy or n-mixture?
occ <- 1
if (occ == 1) {
  y[y > 0] <- 1
} else {
  y <- y
}

data <- list(
  N = nspec, # number of species (integer)
  J = J, # number of sites (integer)
  K = K, # number of visits to each site j (vector [1:80])
  S = nyear, # number of years
  y = y, # 4D array of detection/nondetection data with structure [1:J,1:K,1:N] (in other words, detection matrices by site [J KxN matrices], visit [K JxN matrices], and species [N JxK matrices]) y
  JDay = date1, # a JxK matrix for standardized/centered JDay values for each visit k to site j with no NA present
  JDay2 = date2,
  Time = time1, # a JxK matrix for standardized/centered Time values for each visit k to site j
  Time2 = time2, # same as above with quadratic term for Time
  # Cover=as.vector(Cover), # a vector of length J with standardized/centered values for % cover (1ha)
  # Height=as.vector(Height),
  # Height2=as.vector(Height2)# a vector of length J with standardized/centered values for mean height (1ha)
  # NBR=as.vector(NBR) # a vector of length J with standardized/centered values for mean NBR (1ha)
  El = as.vector(El) # a vector of length J with standardized/centered values for mean elevation (1ha)
  # and any other variables we'd want in either part of the model (aspect, etc... but just starting with these for now)
)


# JAGS code ---------------------------------------------------------------

cat("
    model{

    # COMMUNITY hyperpriors
    omega ~ dunif(0,1) # best practices suggested by Guillera-Arroita et al 2018 suggest dbeta(0.001,0)) but this doesnt work for me

    u.mean ~ dunif(0,1) # Mike Meredith whoever he is prefers dbeta(1,1) to dunif(0,1) for uniform priors
    mu.u <- log(u.mean) - log(1-u.mean) # no idea why this is the way it is
    tau.u ~ dgamma(0.1, 0.1) # Mike Meredith also hates gamma distributions for tau...

    v.mean ~ dunif(0,1)
    mu.v <- log(v.mean) - log(1-v.mean)
    tau.v ~ dgamma(0.1, 0.1)

    # community-level priors for each of the explanatory variables
    mu.aD ~ dnorm(0, 0.001)
    mu.aD2 ~dnorm(0,0.001)
    mu.aT ~ dnorm(0,0.001)
    mu.aT2 ~ dnorm(0,0.001)
    mu.bEl ~ dnorm(0,0.001)

    tau.aD ~ dgamma(0.1,0.1)
    tau.aD2 ~ dgamma(0.1,0.1)
    tau.aT ~ dgamma(0.1,0.1)
    tau.aT2 ~ dgamma(0.1,0.1)
    tau.bEl ~ dgamma(0.1,0.1)

    w[i] ~ dbern(omega) # community-level hyperprior (draw from omega)
    u[i] ~ dnorm(mu.u, tau.u) # occurrence process intercept
    v[i] ~ dnorm(mu.v, tau.v) # detection process intercept

    # DETECTION process
    aDay[i] ~ dnorm(mu.aD, tau.aD)  # effect of date on detection
    aTime[i] ~ dnorm(mu.aT, tau.aT)  # effect of time on detection
    aDay2[i] ~ dnorm(mu.aD2, tau.aD2)  # squared effect of survey date on detection
    aTime2[i] ~ dnorm(mu.aT2, tau.aT2)  # squared effect of survey time on detection
    # OCCURRENCE process
    bEl[i] ~ dnorm(mu.bEl, tau.bEl)  # effect of elevation on occupancy

    for (i in 1:N) {   # loop over species (N)

    for (j in 1:J){    #loop over j sites (J)

    for (s in 1:S) {   #loop over s years (S)

    z[j,i] ~ dbern(mu.psi[j,i]) #state model
    logit(psi[j,i]) <- u[i] + bEl[i]*El[j]
    mu.psi[j,i] <- psi[j,i] * w[i]

    for (k in 1:K[j]) {  #loop over k visits at each point j
    y[j,k,s,i] ~ dbern(mu.p[j,k,s,i]) # detection model
    mu.p[j,k,s,i] <- p[j,k,s,i] * z[j,s,i]
    logit(p[j,k,s,i]) <- v[i] + aDay[i]*JDay[j,k,s] + aTime[i]*Time[j,k,s]

    } # end of year loop
    } # end of visit loop
    } # end of site loop
    } # end of species loop

    # DERIVED quantities
    for (i in 1:N) {  # for each species i, number of occupied sites
      Nocc[i] <- sum(z[,i])
    }
    for (j in 1:J) {  # at each site j, number of species present
      Nsite[j] <- sum(z[j,])
    }


    } # end of model loop
    ", file = "models/jags/Community_Occupancy_trial_DateTime_El.txt")


params <- c(
  "omega", "u", "v",
  "aDay", "aDay2", "aTime", "aTime2",
  "bEl",
  "Nocc", "Nsite"
)

# initial values

zst <- array(1, dim = c(J, S, N)) # for now, starting everything as present (following K+R Ch 11.7.2)
wst <- rep(1, N)
omegaGuess <- runif(1, 0, 1) # not sure where to set this one

inits <- function() {
  list(
    omega = omegaGuess,
    z = zst,
    w = wst,
    u = rnorm(N),
    v = rnorm(N),
    bEl = rnorm(N)
  )
} # these are initial values for jags to start looking for estimates

# run the model

occujags <- jags(data = data, model.file = "models/jags/Community_Occupancy_trial_DateTime_El.txt", inits = inits, parameters.to.save = params, n.chains = 3, n.iter = 10000, n.burnin = 1000)
