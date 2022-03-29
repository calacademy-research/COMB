library(readr)
library(spBayes)
library(classInt)
library(RColorBrewer)
library(MBA)
library(fields)
library(geoR)


year_pt_counts <- read_csv("~/Documents/COMB_DATA/bq-results-20220323-204356-1648068246476.csv")
View(year_pt_counts)

hamfly_df <- subset(subset(year_pt_counts, (bird == "hamfly") & (year == 2020)), !is.na(latitude))
hamfly_ptct <- hamfly_df$detection_count
hamfly_maxlgt <- hamfly_df$max_logit

coords <- as.matrix(hamfly_df[, c("longitude", "latitude")])

x.res <- 100
y.res <- 100
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))

par(mfrow=c(2,1))
## Point Count surface plot
surf1 <- mba.surf(cbind(coords, hamfly_ptct), no.X = x.res, no.Y = y.res, 
                 h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf1, xaxs = "r", yaxs = "r", xlab = "Lon", ylab = "Lat", 
           col = col.br(25), main = "Point Count - Hammonds Flycatcher")
contour(surf, add=TRUE)

## Max logit surface plot
surf2 <- mba.surf(cbind(coords, hamfly_maxlgt), no.X = x.res, no.Y = y.res, 
                 h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf2, xaxs = "r", yaxs = "r", xlab = "Lon", ylab = "Lat",
           col = col.br(25), main = "Max Logit - Hammonds Flycatcher")
contour(surf2, add=TRUE)


#coords <- as.matrix(dat[,c("X","Y")]) 
#nut.names <- c("Ca","K","Mg") 

cor(hamfly_maxlgt,hamfly_ptct)

max <- 0.2*max(as.matrix(dist(hamfly_df[,c("longitude", "latitude")])))

vario_ptct <- variog(coords=coords,data=log(hamfly_ptct+1), 
                  uvec=(seq(0,max, length=24))) 

ml <- likfit(coords=coords,data=log(hamfly_ptct+1), ini = c(1, 0.5))
summary(ml)
log_ptct_hamfly <- log(hamfly_ptct+1)
  
fit <-variofit(vario_ptct, 
               ini.cov.pars=c(5,0.01), 
               cov.model="exponential", 
               minimisation.function="nls", 
               weights="equal",
               fix.nugget = TRUE)
plot(vario, pch=19, main="Max Logit Hammonds Flycatcher") 
lines(fit)
abline(h=fit$nugget, col="blue")
abline(h=fit$cov.pars[1]+fit$nugget, col="green") 
abline(v=-log(0.05)*fit$cov.pars[2], col="red3")

df <- data.frame(log_ptct_hamfly,hamfly_maxlgt)  
q <-2 
n.samples <- 10000
n.ltr <- q*(q+1)/2

# Multivariate model for the joint distribution of point count and max logit
ptct.maxlog.spMvLM <- spMvLM(list(log_ptct_hamfly~1,hamfly_maxlgt~1),
                             coords=coords,
                             data=df,
                             starting= list("beta"=rep(1,q),
                                            "phi" = rep(0.01,q),
                                            "A" = rep(0.1,n.ltr),
                                            "Psi" = rep(0.1,q)),
                                            tuning = list("phi" = rep(0.1,q),
                                                          "A" = rep(0.001,n.ltr),
                                                          "Psi" = rep(0.1,q)),
                                            priors = list("phi.Unif" = list(rep(3/1000,q),
                                                                            rep(3/100,q)),
                                                          "K.IW" = list(q+1, diag(0.1,q)),
                                                          "Psi.IG" = list(rep(2,q),rep(0.1,q))),
                                            cov.model = "exponential",
                                            n.samples = n.samples,
                                            verbose=TRUE,
                                            n.report = 2500
                                            )

burn.in <- floor(0.75*n.samples)
ptct.maxlog.spMvLM <- spRecover(ptct.maxlog.spMvLM, start = burn.in, verbose = FALSE)
summary(ptct.maxlog.spMvLM$p.theta.recover.samples)

## Simple univariate regression model with spatially-varying coefficients
priors <- list("phi.Unif"=list(rep(3/1000,q),
                               rep(3/100,q)),
               "K.IW"=list(q, diag(rep(1,q))),
               "tau.sq.IG"=c(2, 1))
starting <- list("phi"=rep(0.01,q), "A" = rep(0.1,n.ltr), "tau.sq"=1)
tuning <- list("phi"=rep(0.1,q), "A"=rep(0.001, n.ltr), "tau.sq"=0.01)
spatiallm1 <- spSVC(hamfly_maxlgt ~ log_ptct_hamfly, coords=coords, data=df,
               starting=starting, svc.cols=c("(Intercept)","log_ptct_hamfly"),
               tuning=tuning, priors=priors, cov.model="exponential",
               n.samples=10000, n.report=5000)

spatiallm1 <- spRecover(spatiallm1, start=5000, thin=2,  verbose=FALSE)
round(summary(spatiallm1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)


## Simple univariate regression model with non-spatially-varying coefficients
#

# starting values for parameters in spLM()
starting <- list("phi" =3/200, "sigma.sq"=0.08, "tau.sq"=0.02)
# variances for the Metropolis sampler in spLM()
tuning <- list("phi" = 0.1, "sigma.sq"=0.05, "tau.sq"=0.05)
# priors for parameters in spLM(): betas are multivariate Normal
priors.1 <- list("beta.Norm"=list(rep(0,q-1), diag(1000,q-1)),
                  "phi.Unif"=c(3/1500, 3/50), "sigma.sq.IG"=c(2,0.08),
                 "tau.sq.IG"=c(2, 0.02))
# priors for parameters in spLM(): betas are flat
priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                 "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
# function for spatial dependence structure in spLM()
cov.model <- "exponential"
# interval for seeing progress of the sampler in spLM()
n.report <- 5000
# model with first set of priors
spatiallm2 <- spLM(hamfly_maxlgt ~ log_ptct_hamfly - 1, coords=coords, data= df,
                   starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, n.report=n.report)

spatiallm2 <- spRecover(spatiallm2, start=5000, verbose=FALSE)
summary(spatiallm2$p.beta.recover.samples)

hist(log_ptct_hamfly)
hist(hamfly_maxlgt)

spatiallm3 <- spLM(log_ptct_hamfly ~ hamfly_maxlgt - 1, coords=coords, data= df,
                   starting=starting,
                   tuning=tuning, priors=priors.1, cov.model=cov.model,
                   n.samples=n.samples, n.report=n.report)

spatiallm3 <- spRecover(spatiallm3, start=5000, verbose=FALSE)
summary(spatiallm3$p.beta.recover.samples)

