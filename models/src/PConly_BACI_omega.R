library(tidyverse)
library(chron)
library(lubridate)
library(googledrive)
library(here)
library(data.table)
library(reshape2)
library(jagsUI)
library(coda)

source(here("COMB_functions.R"))


# package data ------------------------------------------------------------

PC <- fread(here("point_counts/data_ingest/output/PC_delinted.csv"))
visits <- read_csv("point_counts/data_ingest/output/PC_visit_metadata.csv")

siteList <- unique(visits$point_ID_fk)

# filter to species detected at >3 points
# note: for estimates of richness/community metrics, would be better to include all species, even rare ones. build and run a separate model for richness/diversity estimates.

birdTots <- PC %>%
  group_by(point_ID_fk, birdCode_fk) %>%
  summarise(tot = sum(abun)) %>%
  mutate(detect = if_else(tot > 0, 1, 0)) %>%
  dplyr::select(!tot) %>%
  group_by(birdCode_fk) %>%
  summarise(colSums = sum(detect)) %>% # of points at which the species was detected 
  filter(colSums>3)

dfc <- PC %>% filter(birdCode_fk %in% birdTots$birdCode_fk)

y <- melt(dfc, id.var = c("birdCode_fk", "point_ID_fk", "visit", "year"), measure.var = "abun") %>% acast(point_ID_fk ~ visit ~ year ~ birdCode_fk)

# COVARIATES

dates1 <- visits %>%
  filter(visit < 4) %>%
  mutate(JDay = as.numeric(format(DateTime, "%j"))) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "JDay") %>%
  acast(point_ID_fk ~ visit ~ year)

# standardize dates
mdate <- mean(dates1, na.rm = TRUE)
sddate <- sqrt(var(dates1[1:length(dates1)], na.rm = TRUE))
date1 <- (dates1 - mdate) / sddate

date1[is.na(date1)] <- 0

# TIME matrix

times <- visits %>%
  filter(visit < 4) %>%
  mutate(Time = as.numeric(times(format(DateTime, "%H:%M:%S")))) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "Time") %>%
  acast(point_ID_fk ~ visit ~ year)

# standardize times
mtime <- mean(times, na.rm = TRUE)
sdtime <- sqrt(var(times[1:length(times)], na.rm = TRUE))
time1 <- (times - mtime) / sdtime

time1[is.na(time1)] <- 0

# VEG data

veg <- fread(here("models/input/wide4havars.csv"))

siteData4 <- veg %>%
  dplyr::select(veg_point, avian_point, 
                mean_CanopyBaseHeight_2018_4ha, mean_CanopyBaseHeight_2019_4ha, mean_CanopyBaseHeight_2020_4ha, 
                mean_CanopyHeight_2018_4ha, mean_CanopyHeight_2019_4ha, mean_CanopyHeight_2020_4ha, 
                max_CanopyHeight_2018_4ha, max_CanopyHeight_2019_4ha, max_CanopyHeight_2020_4ha,
                mean_LadderFuelDensity_2018_4ha, mean_LadderFuelDensity_2020_4ha,
                mean_CanopyCover_2018_4ha, mean_CanopyCover_2019_4ha, mean_CanopyCover_2020_4ha, 
                mean_SurfaceFuels_2018_4ha, mean_SurfaceFuels_2019_4ha, mean_SurfaceFuels_2020_4ha, 
                Perc_LTg22mHt_2018_4ha, Perc_LTg22mHt_2020_4ha, Perc_LTg25mHt_2018_4ha, Perc_LTg25mHt_2020_4ha, 
                max_Elevation_NA_4ha, 
                mean_CaplesdNBR_Nov18Nov19_4ha, mean_RAVGcbi4_20182019_4ha) %>%
  arrange(avian_point)

siteData4$RAVG_score <- cut(siteData4$mean_RAVGcbi4_20182019_4ha, breaks = c(-1, 0.5, 2, 3, 4), labels = c(0, 1, 2, 3)) %>% as.factor()

nopc <- setdiff(siteData4$avian_point, siteList) # these three aren't point count locations
setdiff(siteList, siteData4$avian_point) 

siteData4 <- siteData4 %>% filter(avian_point !=0, avian_point !=397, avian_point != 845)

siteList <- sort(unique(dfc$point_ID_fk))

nsite <- length(unique(dfc$point_ID_fk))
nvisit <- length(unique(dfc$visit))
nyear <- length(unique(dfc$year))
nspec <- length(unique(dfc$birdCode_fk))

y <- melt(dfc, id.var = c("birdCode_fk", "point_ID_fk", "visit", "year"), measure.var = "abun") %>% acast(point_ID_fk ~ visit ~ year ~ birdCode_fk)

Cover <- scale(siteData4$mean_CanopyCover_2018_4ha)
#Cover20 <- scale(siteData4$mean_CanopyCover_2020_4ha)
#El <- scale(siteData4$max_Elevation_NA_4ha)
Sev <- round(siteData4$mean_CaplesdNBR_Nov18Nov19_4ha,2)
LgTree <- scale(siteData4$Perc_LTg22mHt_2018_4ha)
#SevS <- (Sev - mean(Sev) / sqrt(var(Sev[1:length(Sev)], na.rm = TRUE)))


BA <- matrix(NA, nsite, nyear)
# colnames(BA) <- 2018:2021
# rownames(BA) <- siteList
BA[, 1:2] <- 0
BA[, 3:4] <- 1

BACI <- BA * as.numeric(round(Sev,2))
BACIb <- BACI
BACIb <- ifelse(BACIb > 0, 1, 0)

spp <- sort(unique(dfc$birdCode_fk)) # VERY important that this is sorted bc this indexes the birds in the later analysis

occ <- 1
if (occ == 1) {
  y[y > 0] <- 1
} else {
  y <- y
}

data <- list(
  nspec = nspec, # number of species (integer)
  nsite = nsite, # number of sites (integer)
  nsurvey = 3, # number of visits to each site
  nyear = nyear, # number of years (integer)
  y = y, # 4D array of detection/nondetection data
  Time = time1, # 3D array of standardized times [nsite x nsurvey x nyear]
  Date = date1,
  LgTree = as.vector(LgTree),
  Cov = as.vector(Cover),
  #Sev = as.vector(Sev),
  #BA = BA,
  BACI = BACIb)


# BACI model specification {JAGS} -----------------------------------------
if (dir.exists(here("models/jags/")) == FALSE) {
  dir.create(here("models/jags"))
}
if (dir.exists(here("models/output/")) == FALSE) {
  dir.create(here("models/output"))
}
# BACI design modified from Cabodevilla et al. 2022

cat("
    model {
    # i = species
    # j = site
    # k = visit
    # t = year

    # OCCURRENCE PRIORS
      # community hyperpriors
    b0.mean.b ~ dbeta(1, 1)
    b0.mean.a ~ dbeta(1, 1)
    mu.b0.b <- logit(b0.mean.b) # separate intercepts before & after fire
    mu.b0.a <- logit(b0.mean.a)
    sd.b0 ~ dunif(0, 5)
    tau.b0 <- 1/sd.b0^2

    mu.bHt ~ dunif(-5, 5)
    sd.bHt ~ dunif(0, 5)
    tau.bHt <- 1/sd.bHt^2

    mu.bCov ~ dunif(-5, 5)
    sd.bCov ~ dunif(0, 5)
    tau.bCov <- 1/sd.bCov^2

    mu.bYr ~ dunif(-5, 5)
    sd.bYr ~ dunif(0, 5)
    tau.bYr <- 1/sd.bYr^2
    
    a0.mean ~ dbeta(1, 1)
    mu.a0 <- logit(a0.mean) 
    sd.a0 ~ dunif(0, 5)
    tau.a0 <- 1/sd.a0^2

    mu.aD ~ dunif(-5, 5)
    sd.aD ~ dunif(0, 5)
    tau.aD <- 1/sd.aD^2

    mu.aT ~ dunif(-5, 5)
    sd.aT ~ dunif(0, 5)
    tau.aT <- 1/sd.aT^2

  # species-level priors

    for (i in 1:nspec) {

        b0[1,i] ~ dnorm(mu.b0.b, tau.b0)
        b0[2,i] ~ dnorm(mu.b0.a, tau.b0)
        bHt[i] ~ dnorm(mu.bHt, tau.bHt)
        bCov[i] ~ dnorm(mu.bCov, tau.bCov)
        bYr[i] ~ dnorm(mu.bYr, tau.bYr)

        a0[i] ~ dnorm(mu.a0, tau.a0)
        aTime[i] ~ dnorm(mu.aT, tau.aT)
        aDate[i] ~ dnorm(mu.aD, tau.aD)

      for (j in 1:nsite) {

        for (t in 1:nyear) {

          logit(psi[j,t,i]) <- b0[BACI[j,t]+1,i] + bHt[i]*LgTree[j] + bCov[i]*Cov[j] + bYr[i]*t

          z[j,t,i] ~ dbern(psi[j,t,i])

          for (k in 1:nsurvey) {

            mu.y[j,k,t,i] <- p[j,k,t,i] * z[j,t,i]
            y[j,k,t,i] ~ dbern(mu.y[j,k,t,i])
            logit(p[j,k,t,i]) <- a0[i] + aTime[i]*Time[j,k,t] + aDate[i]*Date[j,k,t]

                } #k
            } #t
        } #j
    } #i

    # Derived Quantities
      # effect of before/after
      for (i in 1:nspec) {
      effect.ba.sp[i] <- b0[2,i] - b0[1,i]
        for (t in 1:nyear) {
        Nocc[i,t] <- sum(z[,t,i]) # number of points occupied by each species in each yr
        }
      }
      
      for (j in 1:nsite) {
        for (t in 1:nyear) {
          Nsite[j,t] <- sum(z[j,t,]) # species richness (# spp) at each point in each yr
        }
      }
      
    } # end of model loop
    ", file = "models/jags/BACI_occupancy_omega.txt")

# set model params + run ------------------------------------------------------------
params <- c( # I have to run z, psi, and p individually in separate runs or my computer crashes
  # "a0", "aTime", "aDate",
  "b0", "bYr", "bCov", "bHt",
  "effect.ba.sp", "Nsite", "Nocc"
  #"psi", "p", "z"
)

# initial values
zst <- apply(data$y, c(1, 3, 4), max, na.rm = TRUE)

inits <- function() {
  list(z = zst)
} # these are initial values for jags to start looking for estimates

BACI_Out <- jags(data = data, model.file = "models/jags/BACI_occupancy_omega.txt", inits = inits, parameters.to.save = params, n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 3)
out1 <- out
out <- as.data.frame(BACI_Out$summary[1:775, ]) # check length of summary before saving
out$overlap0 <- as.factor(out$overlap0)
#write_csv(out, "models/output/BACI_out.csv")

#save(BACI_Out, out, data, BACI, siteList, Ht, spp, file = "models/output/BACIdata.RData") # need to git ignore this


# output exploration ------------------------------------------------------

library(wesanderson)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(sjPlot)

# Year effect
year <- out[91:135, ]
year$spp <- spp

ggplot(year) +
  geom_point(aes(x = reorder(spp, -mean), y = mean, color = overlap0)) +
  geom_errorbar(aes(x = reorder(spp, -mean), ymin = `2.5%`, ymax = `97.5%`, color = overlap0)) +
  geom_linerange(size = 2, aes(x = reorder(spp, -mean), ymin = `25%`, ymax = `75%`, color = overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y = "posterior estimate of year effect", x = "species code", title = "Effect of survey year on occurrence") +
  theme_classic() +
  scale_color_manual(values = wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50, hjust = 1))

# Canopy Association
cover <- out[136:180, ]
cover$spp <- spp
ggplot(cover) +
  geom_point(aes(x = reorder(spp, -mean), y = mean, color = overlap0)) +
  geom_errorbar(aes(x = reorder(spp, -mean), ymin = `2.5%`, ymax = `97.5%`, color = overlap0)) +
  geom_linerange(size = 2, aes(x = reorder(spp, -mean), ymin = `25%`, ymax = `75%`, color = overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y = "posterior estimate of Canopy Cover effect", x = "species code", title = "Association of canopy cover and bird occurrence") +
  theme_classic() +
  scale_color_manual(values = wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50, hjust = 1))

# Large Tree association

lgtrees <- out[181:225, ]
lgtrees$spp <- spp
ggplot(lgtrees) +
  geom_point(aes(x = reorder(spp, -mean), y = mean, color = overlap0)) +
  geom_errorbar(aes(x = reorder(spp, -mean), ymin = `2.5%`, ymax = `97.5%`, color = overlap0)) +
  geom_linerange(size = 2, aes(x = reorder(spp, -mean), ymin = `25%`, ymax = `75%`, color = overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y = "posterior estimate of Large Tree covariate", x = "species code", title = "Association of large trees and bird occurrence") +
  theme_classic() +
  scale_color_manual(values = wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50, hjust = 1))

# Fire effect (b0[2]-b0[1])
guilds <- read_csv("models/input/CaplesGuilds.csv")
guilds$birdCode_fk[guilds$birdCode_fk=="AUWA"] <- "YRWA" # hacky fix for now, will fix csv later...
guilds$birdCode_fk[guilds$birdCode_fk=="ORJU"] <- "DEJU"
guilds <- guilds %>% filter(birdCode_fk %in% spp) %>% arrange(birdCode_fk)

fire.fx <- out[226:270,]    
fire.fx$spp <- spp
fire.fx <- fire.fx %>% left_join(guilds, by=c("spp"="birdCode_fk")) 
for_cols <- data.frame(fora_sub_1 = (factor(c("Air", "Bark", "Foliage", "Ground"), levels = c("Air", "Bark", "Foliage", "Ground"))),
                       hue = c("lightblue", "yellow", "darkgreen","brown"))
nest_cols <- data.frame(nest_loc_1 = factor(c("Woody Upper Canopy", "Woody Lower Canopy", "Snag", "Shrub", "Ground"), levels = c("Woody Upper Canopy", "Woody Lower Canopy", "Snag", "Shrub", "Ground")),
                        hue = c("darkgreen", "green", "black","orange", "brown"))

f2 <- fire.fx %>% filter(nest_loc_1 !="Cliff")

firenest <- ggplot(f2) + 
  geom_rect(data = nest_cols, aes(fill = hue),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_point(aes(x = reorder(spp, -mean), y = mean)) +
  geom_errorbar(aes(x = reorder(spp, -mean), ymin = `2.5%`, ymax = `97.5%`, linetype=overlap0, width=0.4)) +
  geom_linerange(size = 2, aes(x = reorder(spp, -mean), ymin = `25%`, ymax = `75%`)) +
  geom_hline(yintercept = 0) +
  labs(y = "posterior estimate of fire effect", x = "species code") +
  scale_color_manual(values = c("black", "#3D3D3D")) +
  scale_y_continuous(limits=c(-4,7), breaks = seq(-7,7,1)) +
  scale_linetype_manual(values=c("0"=1,"1"=2)) +
  facet_grid(rows = vars(nest_loc_1), scales = "free_y", space= "free_y") +
  theme(legend.position = "none", strip.text.y = element_text(size = 8, family="Helvetica"), text = element_text(size=10, family = "Helvetica")) +
  coord_flip()

firefor <- ggplot(fire.fx) + 
  geom_rect(data = for_cols, aes(fill = hue),
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_errorbar(aes(x = reorder(spp, -mean), ymin = `2.5%`, ymax = `97.5%`, linetype=overlap0)) +
  geom_linerange(size = 2, aes(x = reorder(spp, -mean), ymin = `25%`, ymax = `75%`)) +
  geom_point(aes(x = reorder(spp, -mean), y = mean)) +
  geom_hline(yintercept = 0) +
  labs(y = "posterior estimate of fire effect", x = "species code") +
  scale_color_manual(values = c("black", "#3D3D3D")) +
  scale_fill_identity() +
  scale_linetype_manual(values=c("0"=1,"1"=2)) +
  theme(legend.position = "none") +
  scale_y_continuous(limits=c(-4,7), breaks = seq(-4,7,1)) +
  facet_grid(rows = vars(fora_sub_1), scales = "free_y", space= "free_y") +
  coord_flip() +
  theme(strip.text.y = element_text(size = 10, face = "bold", family="Helvetica"), text = element_text(size=10, family = "Helvetica"))

## out[271:315, ] # Nsite

# Richness estimates ------------------------------------------------------
# plot of mean sprich (Nsite) for burned vs unburned sites (by burn severity class?
rich <- out[271:594, ]
rich$year <- as.vector(sapply(c("2018", "2019", "2020", "2021"), function(x) rep(x, 81)))
rich$site <- rep(siteList, 4)
rich$RAVG <- rep(as.factor(siteData4$RAVG_score),4)
rich$TallTree <- rep(siteData4$Perc_LTg22mHt_2018_4ha, 4)
rich$TreeBin <- cut(rich$TallTree, breaks = c(-1, 0.02, 0.1, 0.2, 0.55))
rich$Burned <- ifelse(rich$RAVG=="0", "0", "1")
rich$BA <- ifelse(rich$year==2018 | rich$year==2019, "0", "1")

groupmeans <- rich %>% group_by(year, RAVG) %>% summarize(groupmean = mean(mean), groupse = sd(mean))
BAmeans <- rich %>% group_by(BA, RAVG) %>% summarise(BAmean = mean(mean))

# by severity only 
sevCol <- c("darkgreen", "#FFAA00", "#FF6600", "#FF0000")

# retain individual points
ggplot() +
  geom_point(data=rich, aes(x = year, y = mean, group = site, color = RAVG), alpha=.6) +
  geom_line(data=rich, aes(x = year, y = mean, group = site, color = RAVG), alpha=.6) +
  geom_linerange(data=rich, aes(x = year, ymin = `2.5%`, ymax = `97.5%`, color = RAVG), alpha=.6) +
  geom_line(data=groupmeans, aes(x=year, y=groupmean, group = NA), color="black", size=1.2) +
  geom_linerange(data=groupmeans, aes(x = year, ymin = groupmean-groupse, ymax = groupmean+groupse), color="black", size=1.2) +
  #geom_linerange(data=rich, aes(x = year, ymin = `25%`, ymax = `75%`, color = RAVG)) +
  facet_wrap(~RAVG) +
  labs(y = "posterior estimate of plot-level species richness", x = "year", title = "per-point richness over time relative to burn severity") +
  theme_classic() +
  scale_color_manual(values=sevCol) +
  scale_fill_manual(values=sevCol)

ggplot(rich) +
  geom_point(aes(x = year, y = mean, group = site, color = RAVG), alpha=.7) +
  geom_line(aes(x = year, y = mean, group = site, color = RAVG), alpha=.7) +
  geom_linerange(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, color = RAVG), alpha=.7) +
  geom_boxplot(aes(x=year, y=mean), color="black", alpha=0.4) +
  #geom_linerange(data=rich, aes(x = year, ymin = `25%`, ymax = `75%`, color = RAVG)) +
  facet_wrap(~RAVG) +
  labs(y = "posterior estimate of plot-level species richness", x = "year", title = "") +
  theme_classic() +
  scale_color_manual(values=sevCol) +
  scale_fill_manual(values=sevCol)

# by year 
rich %>% group_by(year) %>% summarize(groupmean = mean(mean), groupse = sd(mean))

ggplot(rich, aes(x = mean, y = year, fill = BA)) +
  geom_density_ridges(alpha=0.5, stat="binline", bins=80) +
  theme_classic() + 
  theme(legend.position = "none") +
  labs(x = "posterior mean species richness", y="year") +
  scale_fill_manual(values=sevCol)

ggplot(rich, aes(x = mean, y = year, fill = BA)) +
  geom_density_ridges(alpha=0.5) +
  geom_segment(aes(x = 17.6 , y = 1, xend = 17.6, yend = 4.8), color=sevCol[1], size=1) +
  geom_segment(aes(x = 16.8 , y = 2, xend = 16.8, yend = 5.0), color=sevCol[1], size=1) +
  geom_segment(aes(x = 19.9 , y = 3, xend = 19.9, yend = 5.2), color=sevCol[2], size=1) +
  geom_segment(aes(x = 20.0 , y = 4, xend = 20.0, yend = 5.5), color=sevCol[2], size=1) +
  theme_classic() + 
  theme(legend.position = "none") +
  labs(x = "distribution of posterior mean species richness", y="year") +
  scale_fill_manual(values=sevCol)

# by tree size and severity
bin.labs <- c("0-2%", "2-10%", "11-20%", "20-53%")
names(bin.labs) <- c("(-1,0.02]", "(0.02,0.1]", "(0.1,0.2]", "(0.2,0.55]")

ggplot(rich) +
  geom_line(aes(x = year, y = mean, group = site, color = RAVG)) +
  facet_wrap(~TreeBin, labeller = labeller(TreeBin = bin.labs)) +
  labs(y = "posterior estimate of plot-level species richness", x = "year", title = "Species Richness over time relative to large trees and burn severity") +
  theme_classic() +
  scale_color_manual(values=sevCol)

# Number of sites occupied per yr by each bird ----------------------------

nsitesocc <- out[595:774,]
nsitesocc$spp <- rep(spp,4)
nsitesocc$year <- rep(c("2018", "2019", "2020", "2021"), each=45)

summary(nsitesocc$mean)
hist(nsitesocc$mean)

commonbins <- nsitesocc %>% group_by(spp) %>% 
  summarise(meansites=mean(mean)) %>%
  mutate(common = as.factor(cut(meansites, seq(0,85,5))))

nsitesocc <- left_join(nsitesocc, commonbins)

nsitesocc$BA <- as.factor(ifelse(nsitesocc$year>2019, 1,0))

ggplot(nsitesocc) +
  geom_point(aes(x=year, y=mean, color=spp)) +
  geom_line(aes(x=year, y=mean, group=spp, color=spp)) +
  geom_errorbar(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, color=spp)) +
  facet_wrap(~reorder(spp, -mean)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
