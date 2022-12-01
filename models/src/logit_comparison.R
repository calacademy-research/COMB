#logit_comparison.R
#graphic to show logit distributions
library(ggplot2)
library(dplyr)
library(ggarrange)
library(ggpubr)
library(data.table)

dataML_tall <- fread(paste0(here, "/acoustic/data_ingest/output/dataML_tall.csv"))

focsp1 <- "NOFL"

dataML_tall %>%
  # filter(Date_Time =="2020-06-12 05:30:00") %>%
  # filter(point == "408") %>%
  dplyr::arrange(Start_Time) %>%
  dplyr::filter(logit < 5) %>%
  mutate(focal_species = if_else(species == focsp1, "yes", "no")) %>% 
  ggplot(., aes(x=logit, color = focal_species)) +
  geom_density(adjust = .25) +
  # ylim(1, 200) +
  scale_y_continuous(trans = "log") +
  labs(title=paste("Logits for a single species",
                   focsp1, "versus all others")) +
  # Add mean line
  geom_vline(data = . %>%
               group_by(focal_species) %>%
               summarize(focal_species_mean = mean(logit)),
             mapping = aes(xintercept = focal_species_mean),
             color = c(2,3), linetype="dashed", size=1)

#a function to compare two species logits
logitPlot <- function(dataML_tall, focsp1, focsp2){
  # focsp1 <- "AMDI" 
  # focsp2 <- "YRWA"
  #get dataML (see appropriate 'output' ... )
  dataML_tall %>%
    #filter(Date_Time =="2020-06-12 05:30:00") %>%
    #filter(point == "408") %>%
    arrange(Start_Time) %>%
    mutate(focal_species = if_else(species == focsp1, "yes", "no")) %>%
    ggplot(., aes(x=logit, color = focal_species)) +
    geom_density(adjust = .5) +
    # xlim(-2, 4) +
    labs(title=paste("Logits for a single species",
                     focsp1, "versus all others")) +
    # Add mean line
    geom_vline(data = . %>%
                 group_by(focal_species) %>%
                 summarize(focal_species_mean = mean(logit)),
               mapping = aes(xintercept = focal_species_mean),
               color = c(2,3), linetype="dashed", size=1) -> pfocsp1
  dataML_tall %>%
    #filter(Date_Time =="2020-06-12 05:30:00") %>%
    #filter(point == "408") %>%
    dplyr::arrange(Start_Time) %>%
    mutate(focal_species = if_else(species == focsp2, "yes", "no")) %>%
    ggplot(., aes(x=logit, color = focal_species)) +
    geom_density(adjust = .5) +
    # xlim(-2, 4) +
    labs(title=paste("Logits for a single species",
                     focsp2, "versus all others")) +
    # Add mean line
    geom_vline(data = . %>%
                 group_by(focal_species) %>%
                 summarize(focal_species_mean = mean(logit)),
               mapping = aes(xintercept = focal_species_mean),
               color = c(2,3), linetype="dashed", size=1) -> pfocsp2
  ggarrange(pfocsp1,pfocsp2)
  # rm(pfocsp1,pfocsp2)
}

logitPlot(dataML_tall,"BBWO","NOFL")
