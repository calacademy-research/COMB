library(ggpubr)
library(tidyverse)


# compare different model "sizes"

load("BBWO2020_allmodels_32files.RData")

latlong <- read.csv("models/input/latlong.csv", header = F, col.names = c("siteid", "lat", "long"))


fullsumm <- data.frame(bbwo20full$summary)
arusumm <- data.frame(bbwo20aru$summary)
pcsumm <- data.frame(bbwo20pc$summary)
noscoresumm <- data.frame(bbwo20noscore$summary)
scoreonlysumm <- data.frame(bbwo20scoreonly$summary)
dayaggsumm <- data.frame(dayagg$summary)
dayaggsumm2 <- data.frame(dayagg2$summary)

fullsumm$model <- "full"
arusumm$model <- "aruonly"
pcsumm$model <- "pconly"
noscoresumm$model <- "noscores"
scoreonlysumm$model <- "scoreonly"
dayaggsumm$model <- "dayagg"
dayaggsumm2$model <- "dayagg3"

fullsumm$measure <- rownames(fullsumm)
arusumm$measure <- rownames(arusumm)
pcsumm$measure <- rownames(pcsumm)
noscoresumm$measure <- rownames(noscoresumm)
scoreonlysumm$measure <- rownames(scoreonlysumm)
dayaggsumm$measure <- rownames(dayaggsumm)
dayaggsumm2$measure <- rownames(dayaggsumm2)

allmods <- rbind(fullsumm, arusumm, pcsumm, noscoresumm, scoreonlysumm) #, dayaggsumm, dayaggsumm2)
allmods$model <- factor(allmods$model, levels=c("full", "pconly", "aruonly", "noscores", "scoreonly"))

# 

# MODEL COMPARISONS -------------------------------------------------------
#(site_pos <- data_frame(naive.pc = apply(data$y.pc[, 1, , ], 1, max, na.rm=T), naive.aru = apply(data$y.aru[,1:aruVisitLimit], 1, max, na.rm=TRUE)))
site_pos <- data_frame(naive.pc.n = rowSums(data$y.pc[,1,,], na.rm=TRUE), naive.aru.n = rowSums(data$y.aru[,], na.rm=TRUE))
site_pos$naive.pc.pa <- ifelse(site_pos$naive.pc.n==0, 0, 1)
site_pos$naive.aru.pa <- ifelse(site_pos$naive.aru.n==0,0,1)

site_pos$both <- site_pos$naive.pc.pa + site_pos$naive.aru.pa
site_pos$either <- ifelse(site_pos$both > 0, 1, 0)

jagsResult_BBWO20$mean$mean_psi # Mean of the posterior of the mean occupancy rate psi across sites from model
site_pos$psi <- jagsResult_BBWO20$mean$psi
site_pos$siteid <- as.numeric(rownames(jagsData$y.ind[,]))
site_pos$species <- speciesCode
site_pos$year <- year
tab <- left_join(latlong, site_pos)
write.csv(tab, file = "models/output/bbwo20_32files.csv")

det_tbl <- data.frame(mean = colMeans(site_pos[3:6], na.rm=T), 
                      model=colnames(site_pos[3:6]), 
                      measure="naivePsi") # how close the pc and aru estimates are to the combined model estimates likely differs by species
# do we have a simulation that gives different z matrices for pc, aru, and show that we return a reasonable psi?
det_tbl
psi_tbl <- psi.all %>% filter(measure=="PropOcc") %>% dplyr::select(mean, model, measure)%>% rbind(det_tbl)

jagsResult$mean$PropOcc

ggplot(site_pos) +
  geom_point(aes(x=naive.aru.pa, y=psi)) +
  coord_flip()

# overall psi -------------------------------------------------------------

psi.all <- allmods %>% filter(measure=="mean_psi" | measure=="PropOcc")
propOccplot <- ggplot(psi.all[psi.all$measure=="PropOcc"]) +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.), linewidth=0.4) +
  geom_hline(yintercept = 0.3827, linetype="dashed") + # naive detections by either source
  geom_hline(yintercept = 0.1234, linetype="twodash", color = "#CD9600") + # naive dets by PC only
  geom_hline(yintercept=0.3456, linetype="twodash", color = "#00BE67") + # naive dets by ARU only
  theme_bw() +
  labs(y="posterior finite-sample occupancy") +
 # facet_wrap(~measure, scales = "free_y", ncol=1) +
  scale_y_continuous(breaks=seq(0.1, 0.9, 0.1))

ggsave(plot = psiplot, filename = "models/figures/psi_allmods_BBWO20_20240220.png",width = 6, height = 5, units = "in", dpi = 300)

postocc.plot <- psi.all %>% filter(measure=="PropOcc") %>%
ggplot() +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.), linewidth=0.4) +
  geom_hline(yintercept = 0.3827, linetype="dashed") + # naive detections by either source
  geom_hline(yintercept = 0.1234, linetype="twodash", color = "#CD9600") + # naive dets by PC only
  geom_hline(yintercept=0.3456, linetype="twodash", color = "#00BE67") + # naive dets by ARU only
  theme_bw() +
  labs(y="posterior finite-sample occupancy") +
  # facet_wrap(~measure, scales = "free_y", ncol=1) +
  scale_y_continuous(breaks=seq(0.1, 0.9, 0.1)) +
  theme(legend.position = "none")

ggsave(plot = postocc.plot, filename = "models/figures/postocc_allmods_BBWO20_20240221.png",width = 4, height = 5, units = "in", dpi = 300)


# effect of fire ----------------------------------------------------------

beta1s <- allmods %>% filter(measure=="beta1")

ggplot(beta1s) +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.)) +
  theme_bw() +
  labs(y="posterior BETA1", title="effect of fire severity on BBWO occurrence 2020")
  

## Result QC checks
print(jagsResult)




colMeans(site_pos[1:5], na.rm=T) # how close the pc and aru estimates are to the combined model estimates likely differs by species
# do we have a simulation that gives different z matrices for pc, aru, and show that we return a reasonable psi?
fire.fxall <- allmods %>% filter(grepl("psi.pred.burn", measure))
fire.fxall$Xburn <- rep(jagsData$Xburn, 5)

supp.labs <- c("Point Count Only", "ARU Only", "PC + ARU detections", "PC + ARU Scores")
names(supp.labs) <- c("pconly", "aruonly", "noscores", "scoreonly")

fire.full.p <- fire.fxall %>% filter(model=="full") %>%
  ggplot() +
  geom_ribbon(aes(x=Xburn, ymin=X2.5., ymax=X97.5.), alpha=0.3) +
  geom_ribbon(aes(x=Xburn, ymin=X25., ymax=X75.), alpha=0.3) +
  geom_line(aes(x=Xburn, y=mean), linewidth=1.6) +
  #facet_wrap(~spp) +
  theme_pubclean() +
  labs(y="predicted occupancy probability (psi)", x = "Fire Severity") +
  theme(legend.position = "none") 

fire.all.p <- fire.fxall %>% filter(model !="full") %>%
ggplot() +
  geom_ribbon(aes(x=Xburn, ymin=X2.5., ymax=X97.5., fill=model), alpha=0.3) +
  geom_ribbon(aes(x=Xburn, ymin=X25., ymax=X75., fill=model), alpha=0.3) +
  geom_line(aes(x=Xburn, y=mean, color=model)) +
  #facet_wrap(~spp) +
  theme_pubclean() +
  labs(y="predicted occupancy probability (psi)", x = "Fire Severity") +
  theme(legend.position = "none") +
  facet_wrap(~model)

firefx.p <- ggarrange(fire.full.p, fire.all.p,
                    #  labels = c("A", "B", "C"),
                      ncol=2,
                      widths = c(1,1)
                    #  common.legend = F
                    )

ggsave(plot = firefx.p, filename = "models/figures/firefx_allmods_BBWO20.png",width = 6, height = 5, units = "in", dpi = 300)


# detection probabilities of each model -----------------------------------

# compare p11 with paru01, paru11 for each model type

detcovars <- allmods %>% filter(measure=="p_aru11" | measure=="p_aru01")

ggplot(detcovars) +
  geom_point(aes(x=model, y=mean, color=measure)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5., color=measure)) +
  theme_bw() +
  labs(y="true (p_aru11) and false (p_aru01) detection rates") +
  scale_color_manual(values=c("darkred", "darkgreen"))



# THRESHOLD EXPT ----------------------------------------------------------

# compare different model "sizes"

load("BBWO2020_threshold_expt.RData")

latlong <- read.csv("models/input/latlong.csv", header = F, col.names = c("siteid", "lat", "long"))


summn2 <- data.frame(bbwo20_n2$summary)
summn1 <- data.frame(bbwo20_n1$summary)
summ0 <- data.frame(bbwo20_0$summary)
summ1 <- data.frame(bbwo20_1$summary)
summ2 <- data.frame(bbwo20_2$summary)

summn2$model <- "neg_2"
summn1$model <- "neg_1"
summ0$model <- "zero"
summ1$model <- "1"
summ2$model <- "2"

summn2$measure <- rownames(summn2)
summn1$measure <- rownames(summn1)
summ0$measure <- rownames(summ0)
summ1$measure <- rownames(summ1)
summ2$measure <- rownames(summ2)

threshmods <- rbind(summn2, summn1, summ0, summ1, summ2)
threshmods$model <- factor(threshmods$model, levels=c("neg_2", "neg_1", "zero", "1", "2"))

# 

# MODEL COMPARISONS -------------------------------------------------------


# detection ---------------------------------------------------------------


detcovars <- threshmods %>% filter(measure=="p_aru11" | measure=="p_aru01" | measure=="mean_psi")

display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 8, name = 'Dark2')


detplot <- ggplot(detcovars) +
  geom_point(aes(x=model, y=mean, color=measure)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5., color=measure)) +
  theme_pubclean() +
  scale_y_continuous(name="true (p_aru11) and false (p_aru01) detection probabilities", sec.axis = sec_axis(trans=~.*1,name="predicted occupancy probability")) +
  labs(x="score threshold used in model") +
  scale_color_brewer(type="qual", palette = 2, direction = -1) +
  theme(legend.position = "bottom",
    axis.title.y = element_text(size=13),
    axis.title.y.right = element_text(color = "#7570B3", size=13)
  )

ggsave(detplot, filename = "models/figures/BBWO20_thresh_det.png",width = 6, height = 5, units = "in", dpi = 300)



# overall psi -------------------------------------------------------------

psi.all <- threshmods %>% filter(measure=="mean_psi")
ggplot(psi.all) +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.), linewidth=0.5) +
  theme_bw() +
  labs(y="posterior p(OCCUPANCY)")


# effect of fire ----------------------------------------------------------

beta1s <- threshmods %>% filter(measure=="beta1")

ggplot(beta1s) +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.)) +
  theme_bw() +
  labs(y="posterior BETA1", title="effect of fire severity on BBWO occurrence 2020")


fire.fxall.t <- threshmods %>% filter(grepl("psi.pred.burn", measure))
fire.fxall.t$Xburn <- rep(jagsData$Xburn, 5)



fire.fxall.t  %>%
  ggplot() +
  geom_ribbon(aes(x=Xburn, ymin=X2.5., ymax=X97.5., fill=model), alpha=0.3) +
  geom_ribbon(aes(x=Xburn, ymin=X25., ymax=X75., fill=model), alpha=0.3) +
  geom_line(aes(x=Xburn, y=mean, color=model)) +
  #facet_wrap(~spp) +
  theme_pubclean() +
  labs(y="predicted occupancy probability (psi)", x = "Fire Severity") +
  theme(legend.position = "none") +
  facet_wrap(~model)


agg <- read_csv("models/input/file_logit_agg.csv")
bbwoscores <- agg %>% filter(species=="BBWO", year(Date_Time)==2020)

test <- bbwoscores %>% group_by(point) %>% summarise(mean_max=mean(max_logit)) %>% left_join(site_pos, by=c("point"="siteid")) %>% arrange(psi) #%>% View()

ggplot(test) +
  geom_histogram(aes(mean_max, group=as.factor(both), fill=as.factor(both)), alpha=0.7) 



# Pileated Woodpecker -----------------------------------------------------

piwosumm$model <- "piwo_full"

piwosumm$measure <- rownames(piwosumm)

fire.fx.piwo <- piwosumm %>% filter(grepl("psi.pred.burn", measure))
fire.fx.piwo$Xburn <- rep(jagsData$Xburn, 1)



fire.fx.piwo %>%
  ggplot() +
  geom_ribbon(aes(x=Xburn, ymin=X2.5., ymax=X97.5., fill=model), alpha=0.3) +
  geom_ribbon(aes(x=Xburn, ymin=X25., ymax=X75., fill=model), alpha=0.3) +
  geom_line(aes(x=Xburn, y=mean, color=model)) +
  #facet_wrap(~spp) +
  theme_pubclean() +
  labs(y="predicted occupancy probability (psi)", x = "Fire Severity") +
  theme(legend.position = "none")


# Common Nighthawk --------------------------------------------------------

conisumm <- as.data.frame(jagsResult_CONI0$summary)
conisumm$model <- "coni_full"

conisumm$measure <- rownames(conisumm)


fire.fx.coni <- conisumm %>% filter(grepl("psi.pred.burn", measure))
fire.fx.coni$Xburn <- rep(jagsData$Xburn, 1)
fire.fx.coni$sevBeta <- 1.94716862

coniplot <- ggplot(fire.fx.coni) +
  geom_ribbon(aes(x=Xburn, ymin=`2.5%`, ymax=`97.5%`), fill="#4A2619", alpha=0.6) +
  geom_ribbon(aes(x=Xburn,  ymin=`25%`, ymax=`75%`), fill="#4A2619", alpha=0.6) +
  geom_line(aes(x=Xburn, y=mean)) +
#  scale_fill_continuous_diverging(palette = "Green-Brown") +
#  scale_color_continuous_diverging(palette = "Green-Brown") +   
  #facet_wrap(~spp) +
  theme_pubclean() +
  labs(y="predicted occupancy probability (psi)", x = "Fire Severity", title="Combined Model") +
  theme(legend.position = "none")

ggsave(plot = coniplot, filename = "models/figures/firefx.coni.png",width = 6, height = 5, units = "in", dpi = 300)
