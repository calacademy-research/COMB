load("models/output/BACI_input.RData")
library(tidyverse)
library(chron)
library(jagsUI)
library(coda)
library(lubridate)
library(googledrive)
library(here)
library(data.table)
library(reshape2)

# BACI design modified from Cabodevilla et al. 2022

cat("
    model {
    # i = species
    # j = site
    # k = visit
    # t = year
    
    # OCCURRENCE PRIORS
    mu.b0.before ~ dnorm(0, 0.5) # separate intercepts for before and after fire
    mu.b0.after ~ dnorm(0, 0.5) 
    tau.b0.spp ~ dgamma(0.1, 0.1)
      
    mu.bHt ~ dnorm(0,0.1)
    tau.bHt ~ dgamma(0.1, 0.1)
    
    mu.bInt ~ dnorm(0,0.1)
    tau.bInt ~ dgamma(0.1, 0.1)
      
    mu.bCov ~ dnorm(0,0.1)
    tau.bCov ~ dgamma(0.1, 0.1)
    
    mu.bYr ~ dnorm(0, 0.1)
    tau.bYr ~ dgamma(0.1, 0.1)
    
    # DETECTION PRIORS
    mu.a0 ~ dnorm(0, 0.1) # intercept for detection model 
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mu.aT ~ dnorm(0,0.1)
    tau.aT ~ dgamma(0.1, 0.1)
      
    mu.aD ~ dnorm(0,0.1)
    tau.aD ~ dgamma(0.1, 0.1)
    
    for (i in 1:nspec) {
    
        b0[1,i] ~ dnorm(mu.b0.before, tau.b0.spp)
        b0[2,i] ~ dnorm(mu.b0.after, tau.b0.spp)
        bHt[i] ~ dnorm(mu.bHt, tau.bHt)
        bCov[i] ~ dnorm(mu.bCov, tau.bCov)
        bInt[i] ~ dnorm(mu.bInt, tau.bInt)
        bYr[i] ~ dnorm(mu.bYr, tau.bYr)

        a0[i] ~ dnorm(mu.a0, tau.a0)
        aTime[i] ~ dnorm(mu.aT, tau.aT)
        aDate[i] ~ dnorm(mu.aD, tau.aD)
    
      for (j in 1:nsite) {
      
        for (t in 1:nyear) {
      
          logit(psi[j,t,i]) <- b0[BACI[j,t]+1,i] + bHt[i]*Ht[j] + bCov[i]*Cov[j] +bInt[i]*BACI[j,4]+1*Ht[j] + bYr[i]*t
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
      }
      # species richness at each site/year
      for (j in 1:nsite) {
        for (t in 1:nyear) {
          Nsite[j,t] <- sum(z[j,t,])
      }
      }
    
    } # end of model loop
    ", file ="models/jags/BACI_occupancy.txt") 

# params + run ------------------------------------------------------------
params <- c( # I have to run z, psi, and p individually in separate runs or my computer crashes
 # "a0", "aTime", "aDate",
  "b0", "bYr", "bCov", "bHt", "bInt",
  "effect.ba.sp", "Nsite"
 # "z", "psi", "p"
  )

# initial values
zst <- apply(y,c(1,3,4),max,na.rm=TRUE)  

inits <- function(){list(z = zst)} #these are initial values for jags to start looking for estimates

BACI_Out <- jags(data = data, model.file = "models/jags/BACI_occupancy.txt", inits=inits, parameters.to.save = params, n.chains = 3, n.iter = 10000, n.burnin = 5000, n.thin = 3) 

out <- as.data.frame(BACI_Out$summary[1:657,]) # check length of summary before saving 
out$overlap0 <- as.factor(out$overlap0)
write_csv(out, "models/output/BACI_out.csv")

save(BACI_Out, out, data, BACI, siteList, Ht, spp, file = "models/output/BACIdata.RData") # need to git ignore this