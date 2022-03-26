
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(chron)
library(jagsUI)
library(coda)
library(lubridate)
library(googledrive)
library(here)
library(wesanderson)
library(reshape2)

rm(list = ls())

# Read in data ------------------------------------------------------------

if (dir.exists(here("models/jags/")) == FALSE) {
  dir.create(here("models/jags"))
}

if (dir.exists(here("models/input/")) == FALSE) {
  dir.create(here("models/input"))
}

if (dir.exists(here("models/output/")) == FALSE) {
  dir.create(here("models/output"))
}

# TODO []: link to Google Drive &/or a link to `read_PC.R`. The Google Drive code in readPC.R does not work when using non-CalAcademy login. For now, doing a hard download.
source(here("COMB_functions.R"))
drive_auth()
# 
drive_sync(here("point_counts/data_ingest/input/"), "https://drive.google.com/drive/u/0/folders/13bXkNCZzFwGC8H4k-CJ4Cf3tYJ3t18zW")
# 
# # Reading inputs -------------------
dfc <- fread(here("point_counts/data_ingest/output/PC_delinted.csv")) # make sure these are the specifications on detection distance, years included, etc. you want for your model and rerun delintPC.R with changes to those filters if necessary

y <- melt(dfc, id.var=c("birdCode_fk", "point_ID_fk", "visit"), measure.var="abun") %>% acast(point_ID_fk ~ visit ~ birdCode_fk)

N <- n_distinct(dfc$birdCode_fk)
spp <- sort(unique(dfc$birdCode_fk))

visits <- read_csv("point_counts/data_ingest/output/PC_visit_metadata.csv")

# Covariates --------------------------------------------------------------

# Visit covariates (for detection model) ----------------------------------

# DATE matrix
dates <- visits %>% dplyr::select(-rich) %>% filter(visit<5) %>%
  mutate(JDay = as.numeric(format(DateTime, "%j"))) %>%
  select(-DateTime) %>%
  group_by(point_ID_fk) %>%
  pivot_wider(names_from = visit, values_from = JDay)
dates1<- as.matrix(dates[,2:length(dates)])

# standardize dates
mdate <- mean(dates1, na.rm=TRUE)
sddate <- sqrt(var(dates1[1:length(dates1)], na.rm=TRUE))
date1 <- (dates1-mdate) / sddate
date2 <- date1*date1
date1 <- as.matrix(date1)
date2 <- as.matrix(date2)

date1[is.na(date1)] <- 0

# TIME matrix

times <- visits %>% dplyr::select(-rich) %>% filter(visit<5) %>%
  mutate(Time = times(format(DateTime,"%H:%M:%S"))) %>%
  select(-DateTime) %>%
  group_by(point_ID_fk) %>%
  pivot_wider(names_from = visit, values_from = Time)
times1 <- as.matrix(times[,2:length(times)])

# standardize times  
mtime <- mean(times1, na.rm=TRUE)
sdtime <- sqrt(var(times1[1:length(times1)], na.rm=TRUE))
time1 <- (times1-mtime) /  sdtime
time2 <- time1*time1
time1 <- as.matrix(time1)
time2 <- as.matrix(time2)

time1[is.na(time1)] <- 0
time2[is.na(time2)] <- 0 # for now, setting these to 0

# vector for number of visits per site j
nvisit <- visits %>% filter(visit<5) %>% group_by(point_ID_fk) %>% summarise(nvisits=max(visit))
K <- as_vector(nvisit[,2])


# Site covariates (occurrence model) --------------------------------------
#TODO [] [DDK]: replace with newly generated spatial data 
# the following two csvs were saved from the 3/15/2021 version of 'Caples Points Metadata.xlsx' from shared GoogleDrive-- in the future could add the googledrive package to read in directly?
env <- read_csv("data/caples_env_1ha_20210315.csv")
utms <- read_csv("data/Wildlife_Sampling_Points.csv")
head(env)

prefire <- env %>% dplyr::select(point_d, mean_Cover_2018_1ha, mean_Height_2018_1ha, mean_NBR_2018_1ha, Perc_NtoE_2020_1ha, Binary_NtoE_2020_1ha, mean_Elevation_2020_1ha) %>% left_join(utms, by=c("point_d"="point_id")) %>% left_join(nvisit, by=c("point_d"="point_ID_fk")) %>% arrange(point_d)

siteList <- unique(dfc$point_ID_fk)

setdiff(prefire$point_d, siteList) # these three aren't point count locations
setdiff(siteList, prefire$point_d) # 1072 is missing in the metadata for some reason 

siteData <- prefire %>% filter(point_d %in% unique(visits$point_ID_fk)) %>% separate(SZ_DNS2, c("Size","Density"), sep=1)

J <- n_distinct(siteData$point_d)

# groundtruth -- looks good
siteData %>% ggplot() + geom_boxplot(aes(x=Size, y=mean_Height_2018_1ha)) + labs(main = "Height and Size", y = "mean tree height (1ha)", x = "Tree Size Class (Becky's data)")
siteData %>% ggplot() + geom_boxplot(aes(x=Density, y=mean_Cover_2018_1ha))+ labs(mean = "Density and Cover", y = "mean cover (1ha)", x = "Tree Density Class (Becky's data)")
siteData %>% ggplot() + geom_boxplot(aes(x=Density, y=mean_NBR_2018_1ha)) + labs(main="Density and NBR", y = "NBR (1 ha)", x = "Tree Density Class (Becky's data)")


Cover <- scale(siteData$mean_Cover_2018_1ha)
Height <- scale(siteData$mean_Height_2018_1ha)
NBR <- scale(siteData$mean_NBR_2018_1ha)
Elev <- scale(siteData$mean_Elevation_2020_1ha)

# occupancy or n-mixture? 
occ=1
if (occ==1) { y[y > 0] <- 1 } else {y=y}

data=list(N=N, # number of species (integer) 
          J=J, # number of sites (integer) 
          K=K, # number of visits to each site j (vector [1:80]) 
          y=y, # 3D array of detection/nondetection data with structure [1:J,1:K,1:N] (in other words, detection matrices by site [J KxN matrices], visit [K JxN matrices], and species [N JxK matrices]) y
          JDay=date1, # a JxK matrix for standardized/centered JDay values for each visit k to site j with no NA present (try means imputation)
          Time=time1, # a JxK matrix for standardized/centered Time values for each visit k to site j (see above re: imputation)
          #Time2=time2, # same as above with quadratic term for Time
          #Cover=as.vector(Cover), # a vector of length J with standardized/centered values for % cover (1ha)
          Height=as.vector(Height),
          #Height2=as.vector(Height2)# a vector of length J with standardized/centered values for mean height (1ha)
          #NBR=as.vector(NBR) # a vector of length J with standardized/centered values for mean NBR (1ha)
          Elev=as.vector(Elev) # a vector of length J with standardized/centered values for mean elevation (1ha)
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
    mu.aT ~ dnorm(0,0.001)
    mu.bHt ~ dnorm(0, 0.001)
    mu.bEl ~ dnorm(0,0.001)


    tau.aD ~ dgamma(0.1,0.1)
    tau.aT ~ dgamma(0.1,0.1)
    tau.bHt ~ dgamma(0.1,0.1) 
    tau.bEl ~ dgamma(0.1,0.1) 


    
    for (i in 1:N) {   # loop over species (N)
    
    w[i] ~ dbern(omega) # community-level hyperprior (draw from omega)
    u[i] ~ dnorm(mu.u, tau.u) # occurrence process intercept
    v[i] ~ dnorm(mu.v, tau.v) # detection process intercept
    
    # DETECTION process
    aDay[i] ~ dnorm(mu.aD, tau.aD)
    aTime[i] ~ dnorm(mu.aT, tau.aT)
    # OCCURRENCE process
    bHt[i] ~ dnorm(mu.bHt, tau.bHt)
    bEl[i] ~ dnorm(mu.bEl, tau.bEl)

    
    for (j in 1:J){   #loop over j sites 
    z[j,i] ~ dbern(mu.psi[j,i]) #state model
    logit(psi[j,i]) <- u[i] + bHt[i]*Height[j] + bEl[i]*Elev[j] 
    mu.psi[j,i] <- psi[j,i] * w[i]   
    
    for (k in 1:K[j]) {  #loop over k visits at each point j
    y[j,k,i] ~ dbern(mu.p[j,k,i]) # detection model 
    mu.p[j,k,i] <- p[j,k,i] * z[j,i]
    logit(p[j,k,i]) <- v[i] + aDay[i]*JDay[j,k] + aTime[i]*Time[j,k]
    
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
    ", file ="models/jags/Community_Occupancy_2019_DateTime_HtElev.txt")


params <- c("omega", "u", "v", 
            "aDay", "aTime",
            "bHt", "bEl",
            "Nocc", "Nsite")        

# initial values

zst <- array(1,dim=c(J,N)) # for now, starting everything as present (following K+R Ch 11.7.2)
wst <- rep(1,N) 
omegaGuess <- runif(1,0,1) # not sure where to set this one

inits <- function(){list(omega = omegaGuess, 
                         z = zst, 
                         w=wst, 
                         u = rnorm(N), 
                         v = rnorm(N), 
                         bHt =rnorm(N), 
                         bEl = rnorm(N))} #these are initial values for jags to start looking for estimates

#code to actually run the model is below

occujags <- jags(data = data, model.file = "models/jags/Community_Occupancy_2019_DateTime_HtElev.txt", inits = inits, parameters.to.save = params, n.chains = 3, n.iter = 10000, n.burnin = 1000) 



# Save model results ------------------------------------------------------

# UNDER CONSTRUCTION: eventually the section below, model exploration, will go in a separate script

plot(occujags) # check traces 
# WHY ARE THESE TWO THINGS DIFFERENT
#apply(plogis(occujags$sims.list$u),2,mean)
#plogis(occujags$summary[2:49,])
# fixed-- bc order of operations is incorrect in the first one
# probability of occupancy for each species
est.occ <- plogis(apply(occujags$sims.list$u,2,mean))
obs.occ <- apply(y,3,sum) # total detections made
# extracting posteriors and confidence intervals
output <- as.data.frame(occujags$summary[1:337,])
output$overlap0 <- as.factor(output$overlap0)

# 
whiskerplot(occujags, parameters = "u")

post_Occ <- output[2:41,] # occupancy estimates
post_Occ$species <- spp
post_Occ$plogis <- plogis(post_Occ$mean)
post_Occ <- post_Occ %>% mutate_at(vars(mean:`97.5%`), plogis)

ggplot(post_Occ) +
  geom_point(aes(x=reorder(species, -plogis), y=plogis, color=overlap0)) +
  geom_errorbar(aes(x=reorder(species,-plogis), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(species,-plogis), ymin=`25%`, ymax=`75%`,color=overlap0)) +
  labs(y="posterior estimate of Occupancy", x = "species code") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1))

output[82:121,] # detection estimates

post_Date <- output[82:121,] # Day on detection for each spp
post_Date$species <- spp
ggplot(post_Date) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=reorder(species, mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(species, mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(species,mean), ymin=`25%`, ymax=`75%`,color=overlap0)) +
  labs(y="posterior estimate of Date on detection", x = "species code") +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none")

## TIME
post_Time <- output[122:161,] # Time on detection for each spp
post_Time$species <- spp
ggplot(post_Time) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=reorder(species, mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(species, mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(species,mean), ymin=`25%`, ymax=`75%`,color=overlap0)) +
  labs(y="posterior estimate of Time on detection", x = "species code") +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none")

## HEIGHT 
post_Ht <- output[162:201,] 
post_Ht$species <- spp
ggplot(post_Ht) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=reorder(species, mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(species, mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(species,mean), ymin=`25%`, ymax=`75%`,color=overlap0)) +
  labs(y="posterior estimate of Height", x = "species code") +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none")

## ELEVATION
post_El <- output[202:241,] # Elev "" ""
post_El$species <- spp
ggplot(post_El) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=reorder(species, mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(species, mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(species,mean), ymin=`25%`, ymax=`75%`,color=overlap0)) +
  labs(y="posterior estimate of Elevation", x = "species code") +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values= wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none")


# [UNDER CONSTRUCTION] more model exploration -----------------------------

# 2022-03-17: need to edit below code for correct indices in occujags

# More common birds have less uncertainty in their occupancy estimates 
# (code based on K+R vol. 1 p 666)
cbind(obs.occu = obs.occ, occujags$summary[338:385,c(1,3,7)])
plot(obs.occ, occujags$summary[338:385,1], xlab="total observed detections", ylab="estimated # sites occupied", frame=F, pch=16)
segments(obs.occ, occujags$summary[338:385,3], obs.occ, occujags$summary[338:385,7])

# Relationship between detection and occupancy probabilities +/- 95% CI

plot(plogis(occujags$summary[2:49,1]), plogis(occujags$summary[50:97,1]), xlab="Occupancy estimate", ylab="Detection estimate", pch=16, xlim=c(-0.01, 1), ylim=c(-0.01, 1), main="As p(detect) increases, certainty in p(occ) increases")

segments(plogis(occujags$summary[2:49,3]), plogis(occujags$summary[50:97,1]), plogis(occujags$summary[2:49,7]), plogis(occujags$summary[50:97,1]), col="grey")
segments(plogis(occujags$summary[2:49,1]), plogis(occujags$summary[50:97,3]), plogis(occujags$summary[2:49,1]), plogis(occujags$summary[50:97,7]), col="grey")

# Caples-level occupancy for each species
whiskerplot(occujags, parameters = "u")
# cover
occujags$summary[146:193,]
whiskerplot(occujags, parameters = "bCov")

# NBR
occujags$summary[194:244,]
whiskerplot(occujags, parameters = "bNBR")

# height 
post_ht <- as.data.frame(occujags$summary[193:250,])
whiskerplot(occujags, parameters = "bHt")
post_ht$species <- sort(unique(PCmat$species))
# Explore derived quantities ----------------------------------------------

# species richness by point
rich = occujags$sims.list$Nsite
mean(rich)
summary(rich)
plot(table(rich))

#Plot mean site richness against one of the covariates to examine how point 
#richness varies as a result of understory foliage
plot(siteData$mean_Height_2018_1ha, apply(rich,2,mean), pch=16, lwd=2, xlab="Mean Tree Height",
     ylab="Point richness", type="p")

plot(siteData$mean_Elevation_2020_1ha, apply(rich,2,mean), pch=16, lwd=2, xlab="Elevation (m)", ylab="Point richness", type="p")

# number of sites occupied by each species
Nocc_spp <- occujags$sims.list$Nocc
mean(Nocc_spp)
summary(Nocc_spp)
plot(table(Nocc_spp))

apply(Nocc_spp,2,mean)

plot(y$birdCode_fk, apply(Nocc_spp,2,mean), pch=16, lwd=2, xlab="Mean Tree Height",
     ylab="Point richness", type="p")
plot(siteData$mean_Elevation_2020_1ha, apply(rich,2,mean), pch=16, lwd=2, xlab="Elevation (m)", ylab="Point richness", type="p")

