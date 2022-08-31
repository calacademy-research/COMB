#graphic to show logit distributions
library(ggplot2)
library(dplyr)
library(ggarrange)
library(ggpubr)

#a function to compare two species logits
logitPlot <- function(dataML, focsp1, focsp2){
pfocsp1 <- NULL 
pfocsp2 <- NULL
#get dataML (see appropriate 'output' ... )
dataML %>%
#  filter(Date_Time =="2020-06-11 09:30:00") %>%
# filter(point == "444") %>%
  arrange(Start_Time) %>%
  mutate(assign(paste0(focsp1)) = if_else(species == focsp1, "yes", "no")) %>%
  ggplot(., aes(x=logit, color = focal_species)) +
  geom_density(adjust = 1) +
  xlim(-2, 4) +
  labs(title=paste("Logits for a single species",
  focsp1, "versus all others")) +
  # Add mean line
  geom_vline(data = . %>%
  group_by(focal_species) %>%
  summarize(focal_species_mean = mean(logit)),
  mapping = aes(xintercept = focal_species_mean),
  color = c(2,3), linetype="dashed", size=1) -> pfocsp1
dataML %>%
  #  filter(Date_Time =="2020-06-11 09:30:00") %>%
  # filter(point == "444") %>%
  arrange(Start_Time) %>%
  mutate(focal_species = if_else(species == focsp2, "yes", "no")) %>%
  ggplot(., aes(x=logit, color = focal_species)) +
  geom_density(adjust = 1) +
  xlim(-2, 4) +
  labs(title=paste("Logits for a single species",
                   focsp2, "versus all others")) +
  # Add mean line
  geom_vline(data = . %>%
             group_by(focal_species) %>%
             summarize(focal_species_mean = mean(logit)),
             mapping = aes(xintercept = focal_species_mean),
             color = c(2,3), linetype="dashed", size=1) -> pfocsp2
ggarrange(pfocsp1,pfocsp2)
rm(pfocsp1,pfocsp2)
}

logitPlot(dataML_top1,"HAWO","PIWO")
