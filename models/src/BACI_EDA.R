load("models/output/BACIdata.RData")

library(wesanderson)
library(tidyverse)
library(AHM)

# Occurrence covariates ---------------------------------------------------

### Year effect -------------------------------------------------------------

out <- BACI_out
out$overlap0 <- as.factor(out$overlap0)

year <- out[97:144,]
year$spp <- spp

ggplot(year) +
  geom_point(aes(x=reorder(spp,-mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(spp,-mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(spp,-mean), ymin=`25%`, ymax=`75%`, color=overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y="posterior estimate of year effect", x = "species code", title="Effect of survey year on occurrence") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1))


### Cover effect ------------------------------------------------------------


cover <- out[145:192,]
cover$spp <- spp
ggplot(cover) +
  geom_point(aes(x=reorder(spp, -mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(spp, -mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(spp,-mean), ymin=`25%`, ymax=`75%`, color=overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y="posterior estimate of Canopy Cover effect", x = "species code", title="Effect of canopy cover on occurrence") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1))


### Large Tree effect -------------------------------------------------------



trees <- out[193:240,]
trees$spp <- spp

# trees %>% filter(spp=="BBWO" | spp=="NOFL" | spp=="MOUQ" | spp=="WETA" |spp=="OSFL" | spp=="WHWO" | spp=="STJA" | spp=="SOGR" | spp=="WEWP" | spp=="HAWO" | spp=="SPTO" | spp=="BHGR" | spp=="DUFL" | spp=="FOSP" | spp=="HEWA" | spp=="AUWA" | spp=="HETH" | spp=="GCKI") %>%
ggplot(trees) +
  geom_point(aes(x=reorder(spp,-mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(spp,-mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(spp,-mean), ymin=`25%`, ymax=`75%`, color=overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y="posterior estimate of Tree Height effect", x = "species code", title="Effect of large trees on occurrence") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1))

tree.spp <- trees %>% filter(overlap0==0)


### Fire effect (b0After - b0Before) ----------------------------------------



fire.fx <- out[289:336,]
fire.fx$spp <- spp

# fire.fx %>% filter(spp %in% tree.spp$spp) %>%
#   mutate(updown = ifelse(mean>0,"up","down")) %>%
ggplot(fire.fx) +
  geom_point(aes(x=reorder(spp,-mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(spp,-mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(spp,-mean), ymin=`25%`, ymax=`75%`, color=overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y="posterior estimate of fire effect", x = "species code", title="Effect of fire on occurrence") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1)) 


### Interaction effect ------------------------------------------------------


int <- out[241:288,]
int$spp <- spp

# fire.fx %>% filter(spp %in% tree.spp$spp) %>%
#   mutate(updown = ifelse(mean>0,"up","down")) %>%
ggplot(int) +
  geom_point(aes(x=reorder(spp,-mean), y=mean, color=overlap0)) +
  geom_errorbar(aes(x=reorder(spp,-mean), ymin=`2.5%`, ymax=`97.5%`, color=overlap0)) +
  geom_linerange(size=2,aes(x=reorder(spp,-mean), ymin=`25%`, ymax=`75%`, color=overlap0)) +
  geom_hline(yintercept = 0) +
  labs(y="posterior estimate of interaction effect", x = "species code", title="Interaction of large trees x fire") +
  theme_classic() +
  scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1)) 

# Trees x Fire Biplot [UNDER CONSTRUCTION] --------------------------------

# colnames(trees) <- paste("tree", colnames(trees), sep = "_")
# colnames(fire.fx) <- paste("fire", colnames(fire.fx), sep = "_")
# 
# d.biplot <- trees %>% 
#   left_join(fire.fx, by=c("tree_spp" = "fire_spp"))
# 
# d.biplot$tree_fx[d.biplot$tree_mean < 0 & d.biplot$tree_overlap0==0] <- "negative"
# d.biplot$tree_fx[d.biplot$tree_mean > 0 & d.biplot$tree_overlap0==0] ="positive" 
# d.biplot$tree_fx[d.biplot$tree_overlap0==1] ="no effect" 
# d.biplot$fire_fx[d.biplot$fire_mean < 0 & d.biplot$fire_overlap0==0] <- "negative"
# d.biplot$fire_fx[d.biplot$fire_mean > 0 & d.biplot$fire_overlap0==0] ="positive" 
# d.biplot$fire_fx[d.biplot$fire_overlap0==1] ="no effect" 
# 
# d.biplot <- d.biplot %>% mutate(sig = case_when(tree_overlap0==0 & fire_overlap0==0 ~ 4,
#                          tree_overlap0==0 & fire_overlap0==1 ~ 3,
#                          tree_overlap0==1 & fire_overlap0==0 ~ 2,
#                          tree_overlap0==1 & fire_overlap0==1 ~ 1))
#          
# head(d.biplot)
# rownames(d.biplot) <- d.biplot$tree_spp
# 
# ggplot(d.biplot, aes(x=tree_mean, y=fire_mean, color=as.factor(sig))) +
#   geom_text(label=rownames(d.biplot), fontface="bold") +
#   #geom_errorbar(aes(x=tree_mean, ymin=`fire_2.5%`, ymax=`fire_97.5%`, color=as.factor(sig))) +
#   #geom_errorbar(aes(y=fire_mean, xmin=`tree_2.5%`, xmax=`tree_97.5%`, color=as.factor(sig), alpha=0.4)) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   scale_colour_manual(values = c("azure4", "darkorange2", "forestgreen", "goldenrod3"),
#                       labels = c("no effect", "effect of fire", "effect of trees", "effect of both"))  +
#   labs(x="effect of large trees on occurrence", y= "effect of fire on occurrence intercept", color="parameter significance:") +
#   theme_classic() +
#   theme(legend.position = "top") -> biplot
# biplot
# 
# bin.labs.f <- c("neg effect of fire", "pos effect of fire")
# bin.labs.t <- c("neg effect of trees", "pos effect of trees")
# names(bin.labs.f) <- c(-1,1)
# names(bin.labs.t) <- c(-1,1)
# 
# biplot +
#   facet_grid(fire_fx ~ tree_fx, labeller = labeller(tree_fx = bin.labs.t, 
#                                                     fire_fx = bin.labs.f)) +
#   theme_bw() +
#   theme(legend.position="top")
# 
# # facet wrap by tree assn
# ggplot(d.biplot) +
#   geom_point(aes(x=reorder(spp,-fire_mean), y=fire_mean, color=fire_overlap0)) +
#   geom_errorbar(aes(x=reorder(spp,-fire_mean), ymin=`fire_2.5%`, ymax=`fire_97.5%`, color=fire_overlap0)) +
#   geom_linerange(size=2,aes(x=reorder(spp,-fire_mean), ymin=`fire_25%`, ymax=`fire_75%`, color=fire_overlap0)) +
#   geom_hline(yintercept = 0) +
#   labs(y="effect of fire on occurrence intercept", x = "species code", title="Effect of fire on occurrence, faceted by tree association") +
#   theme_classic() +
#   scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1)) +
#   facet_grid(tree_fx ~ .)
# 
# ggplot(d.biplot) +
#   geom_point(aes(x=reorder(spp,-tree_mean), y=tree_mean, color=tree_overlap0)) +
#   geom_errorbar(aes(x=reorder(spp,-tree_mean), ymin=`tree_2.5%`, ymax=`tree_97.5%`, color=tree_overlap0)) +
#   geom_linerange(size=2,aes(x=reorder(spp,-tree_mean), ymin=`tree_25%`, ymax=`tree_75%`, color=tree_overlap0)) +
#   geom_hline(yintercept = 0) +
#   labs(y="posterior estimate of tree effect", x = "species code", title="Effect of trees on occurrence, faceted by fire response") +
#   theme_classic() +
#   scale_color_manual(values=wes_palette("BottleRocket1", n = 2)) +
#   theme(legend.position = "none", axis.text.x = element_text(angle = 50,hjust=1)) +
#   facet_grid(fire_fx ~ .)


# Richness estimates ------------------------------------------------------
# plot of mean sprich (Nsite) for burned vs unburned sites (by burn severity class?
rich <- out[337:656,]
rich$year <- as.vector(sapply(c("2018", "2019", "2020", "2021"), function(x) rep(x,80)))
rich$site <- rep(siteList,4)
rich$RAVG <- rep(Sev,4)
rich$TallTree <- rep(siteData4$Perc_LTg22mHt_2018_4ha, 4)
rich$TreeBin <- cut(rich$TallTree, breaks = c(-1, 0.02, 0.1, 0.2, 0.5))

bin.labs <- c("0-2%", "2-10%", "11-20%", "20-50%")
names(bin.labs) <- c("(-1,0.02]", "(0.02,0.1]","(0.1,0.2]", "(0.2,0.5]")

ggplot(rich) +
  geom_line(aes(x=year, y=mean, group=site, color=-RAVG)) +
  facet_wrap(~TreeBin, labeller = labeller(TreeBin = bin.labs)) +
  labs(y="posterior estimate of plot-level species richness", x = "year", title="Species Richness over time relative to large trees and burn severity") +
  theme_classic() +
  scale_color_gradient()

unburned <- rich %>% filter(RAVG=="0") 
ggplot(unburned) +
  geom_line(aes(x=year, y=mean, group=site))

rich %>% filter(RAVG != "0") %>%
ggplot() +
  geom_line(aes(x=year, y=mean, group=site))

rich %>% ggplot() +
  geom_boxplot(aes(x=year,y=mean)) +
  facet_wrap(~RAVG)

# before and after betas
BA <- out[1:96,]
BA$spp <- as.vector(sapply(spp, function(x) rep(x,2)))
BA$prepost <- rep(c("pre", "post"))

pre <- BA %>% filter(prepost=="pre")
post <- BA %>% filter(prepost=="post")
post$firesig <- fire.fx$overlap0
pre$firesig <- fire.fx$overlap0

post2 <- post %>% filter(post$firesig=="0")
pre2 <- pre %>% filter(pre$firesig=="0")

# Posterior Predictions by Species -------------------------------------------------------------
# create new design matrix with "new" covariates using identical scaling as our original data
summary(data$Ht)
o.tree <- seq(min(data$Ht), max(data$Ht),,500)

summary(data$Cov)
o.cov <- seq(min(data$Cov), max(data$Cov),,500)

str(tmp <- BACI_Out$sims.list)
nsamp <- 4998
nsamp

### community mean occupancy predictions ------------------------------------
# 
# predC <- array(NA, dim = c(500,nsamp,6)) 
# for (i in 1:nsamp) {
#   predC[,i,1] <- plogis(tmp$b0[i] + tmp$bHt[i]*o.tree)
#   predC[,i,2] <- plogis(tmp$b0[i] + tmp$bSev[i]*o.sev)
#   predC[,i,3] <- plogis(tmp$b0[i] + tmp$bInt[i]*o.sev*o.tree)
# }
# 
# head(predC)
# pmC <- apply(predC, c(1,3), mean)
# criC <- apply(predC, c(1,3), function(x) quantile(x, prob=c(0.025,0.975)))
# 
# plot(o.tree, pmC[,1], col="blue", lwd=3, type="l", lty=1, frame=F,
#      ylim=c(0,1), xlab="% Large Trees", ylab="Community Mean Occupancy")
# matlines(o.tree, t(criC[,,1]), col="grey", lty=1)
# 
# plot(o.sev, pmC[,2], col="blue", lwd=3, type="l", lty=1, frame=F,
#      ylim=c(0,0.8), xlab="Fire Severity", ylab="Community Mean Occupancy")
# matlines(o.sev, t(criC[,,2]), col="grey", lty=1)
# 
# pmC[,3]
# head(pmC[,3])


### species-specific predictions --------------------------------------------

str(tmp) # this is the sims.list (or just the "naked" mcmc runs)
tmp <- as.array(tmp)
pm <- BACI_Out$mean
cr.lo <- BACI_Out$q2.5
cr.hi <- BACI_Out$q97.5
head(pm)

psi.coef <- cbind(pm$b0[1,], pm$b0[2,], pm$bHt, pm$bCov, pm$bInt, pm$effect.ba.sp) # setting it to b0[2] for the post-fire betas

psi.cilo <- cbind(cr.lo$b0[1,], cr.lo$b0[2,], cr.lo$bHt, cr.lo$bCov, cr.lo$bInt, cr.lo$effect.ba.sp)

psi.cihi <- cbind(cr.hi$b0[1,], cr.hi$b0[2,], cr.hi$bHt, cr.hi$bCov, cr.hi$bInt, cr.hi$effect.ba.sp)

colnames(psi.coef) <- c("bBefore", "bAfter", "bHt", "bCov", "bInt", "bBACI")
colnames(psi.cihi) <- c("bBefore", "bAfter", "bHt", "bCov", "bInt", "bBACI")
colnames(psi.cilo) <- c("bBefore", "bAfter", "bHt", "bCov", "bInt", "bBACI")
rownames(psi.coef) <- spp

simHt <- o.tree
simCov <- o.cov

predS <- array(NA, dim=c(500,48,4))
for(i in 1:48) {
  predS[,i,1] <- plogis(psi.coef[i,1] + psi.coef[i,3]*simHt)
  predS[,i,2] <- plogis(psi.coef[i,1] + psi.coef[i,4]*simCov)
  predS[,i,3] <- plogis(psi.coef[i,1] + 
                          psi.coef[i,3]*simHt + # additive effect of large trees
                          #psi.coef[i,4]*simCov + # effect of cover
                          psi.coef[i,5]*simHt*1) # unburned
  predS[,i,4] <- plogis(psi.coef[i,2] + 
                          psi.coef[i,3]*simHt +
                          #psi.coef[i,4]*simCov +
                          psi.coef[i,5]*simHt*2) # burned
}
# credible intervals [2.5-97.5%]
intS <- array(NA, dim=c(500,48,8))
for(i in 1:48) {
  intS[,i,1] <- plogis(psi.cihi[i,1] + psi.cihi[i,4]*simCov)
  intS[,i,2] <- plogis(psi.cilo[i,1] + psi.cilo[i,4]*simCov)
  intS[,i,3] <- plogis(psi.cihi[i,1] + psi.cihi[i,3]*simHt)
  intS[,i,4] <- plogis(psi.cilo[i,1] + psi.cilo[i,3]*simHt)
  intS[,i,5] <- plogis(psi.cihi[i,2] + psi.cihi[i,5]*simHt*2)
  intS[,i,6] <- plogis(psi.cilo[i,2] + psi.cilo[i,5]*simHt*2)
  intS[,i,5] <- plogis(psi.cihi[i,1] + psi.cihi[i,5]*simHt*1)
  intS[,i,6] <- plogis(psi.cilo[i,1] + psi.cilo[i,5]*simHt*1)
}

plot(simHt, predS[,1,1], lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability")
for (i in 1:48) {
  lines(simHt, predS[,i,1], col=i, lwd=3)
}

# COVER

par(mfrow=c(4,3))
for(i in 1:48) {
  plot(simCov, predS[,i,2], lwd=3, type="l", lty=1, frame=F, ylim = c(0,1),
       xlab="% Canopy Cover", ylab = "pre-fire psi", main=spp[i])
  lines(simCov, intS[,i,1], lty=3)
  lines(simCov, intS[,i,2], lty=3)
}

par(mfrow=c(3,3))
for(i in 1:48) {
  plot(simHt, predS[,i,1], lwd=3, type="l", lty=1, frame=F, ylim = c(0,1), col="gray",
       xlab="% Large Tree Cover", ylab = "occupancy probability", main=spp[i])
  lines(simHt, predS[,i,3], col="forestgreen", lwd=3)
  lines(simHt, predS[,i,4], col="red", lwd=3)
}

# Fox Sparrow
par(mfrow=c(2,3))
plot(simHt, predS[,14,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Fox Sparrow")
lines(simHt, predS[,14,2], col="purple", lwd=3)
lines(simHt, predS[,14,3], col="forestgreen", lwd=3)
lines(simHt, predS[,14,4], col="red", lwd=3)

# GCKI
plot(simHt, predS[,15,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Golden-crowned Kinglet")
lines(simHt, predS[,15,3], col="forestgreen", lwd=3)
lines(simHt, predS[,15,4], col="red", lwd=3)

# GTTO
plot(simHt, predS[,16,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Green-tailed Towhee")
lines(simHt, predS[,16,3], col="forestgreen", lwd=3)
lines(simHt, predS[,16,4], col="red", lwd=3) 

# HETH
plot(simHt, predS[,19,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Hermit Thrush")
lines(simHt, predS[,19,3], col="forestgreen", lwd=3)
lines(simHt, predS[,19,4], col="red", lwd=3) 

# HEWA
plot(simHt, predS[,20,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Hermit Warbler")
lines(simHt, predS[,20,3], col="forestgreen", lwd=3)
lines(simHt, predS[,20,4], col="red", lwd=3) 

# MGWA
plot(simHt, predS[,24,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="MacGillivray's Warbler")
lines(simHt, predS[,24,3], col="forestgreen", lwd=3)
lines(simHt, predS[,24,4], col="red", lwd=3) 

# AUWA
plot(simHt, predS[,3,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Yellow-rumped Warbler")
lines(simHt, predS[,3,3], col="forestgreen", lwd=3)
lines(simHt, predS[,3,4], col="red", lwd=3) 

plot(simHt, predS[,4,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Black-backed Woodpecker")
lines(simHt, predS[,4,3], col="forestgreen", lwd=3)
lines(simHt, predS[,4,4], col="red", lwd=3) 

plot(simHt, predS[,45,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Western Tanager")
lines(simHt, predS[,45,3], col="forestgreen", lwd=3)
lines(simHt, predS[,45,4], col="red", lwd=3) 

plot(simHt, predS[,30,1], col="gray", lwd=3, type="l", lty=1, frame=F, ylim= c(0,1), 
     xlab="% Large Tree Cover", ylab = "Occupancy probability", main="Olive-sided Flycatcher")
lines(simHt, predS[,30,3], col="forestgreen", lwd=3)
lines(simHt, predS[,30,4], col="red", lwd=3) 


# z scores

# psi exploration [UNDER CONSTRUCTION] ---------------------------------------------------------
# not ready to run 
# psi <- BACI_z$summary[15361:30720,]
# p <- BACI_z$summary[30721:76801,]
# z <- BACI_z$summary[1:15360,]
# 
# head(psi)
# 
# for (i in 1:psi.list) {
#   psi.list$i$year <- rep(yearList, 48)
# } 
# 
# # add species names
# psi.list <- split(as.data.frame(psi), rep(1:80))
# names(psi.list) <- siteList
# psi$year <- as.vector(sapply(c("2018", "2019", "2020", "2021"), function(x) rep(x,80)))
# rich$site <- rep(siteList,4)
# # add sites
# # add years