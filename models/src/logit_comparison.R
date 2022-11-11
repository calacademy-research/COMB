#logit_comparison.R
#graphic to show logit distributions
library(ggplot2)
library(dplyr)
library(ggarrange)
library(ggpubr)

dataML_tall <- fread(paste0(here, "/acoustic/data_ingest/output/dataML_tall.csv"))

focsp1 <- "NAWA"

dataML_tall %>%
  # filter(Date_Time =="2020-06-12 05:30:00") %>%
  # filter(point == "408") %>%
  dplyr::arrange(Start_Time) %>%
  dplyr::filter(logit < 5) %>%
  mutate(focal_species = if_else(species == focsp1, "yes", "no")) %>% 
  ggplot(., aes(x=logit, color = focal_species)) +
  geom_density(adjust = .5) +
  # ylim(1, 200) +
  scale_y_continuous(trans = "log") +
  labs(title=paste("Logits for a single species",
                   focsp1, "versus all others")) +
  # Add mean line
  geom_vline(data = . %>%
               group_by(focal_species) %>%
               summarize(focal_species_mean = exp(mean(log(logit+4)))),
             mapping = aes(xintercept = focal_species_mean),
             color = c(2,3), linetype="dashed", size=1)


# model.0 <- lm(log(y - theta.0) ~ x, data=data.df)  

## Not run:  lambda1 <- exp(1); lambda2 <- exp(3)
(phi <- logitlink(-1, inverse = TRUE))
mdata <- data.frame(y1 = rexp(nn <- 1000, lambda1))
mdata <- transform(mdata, y2 = rexp(nn, lambda2))
mdata <- transform(mdata, Y  = ifelse(runif(nn) < phi, y1, y2))
fit <- vglm(Y ~ 1, mix2exp, data = mdata, trace = TRUE)
coef(fit, matrix = TRUE)

# Compare the results with the truth
round(rbind('Estimated' = Coef(fit),
            'Truth' = c(phi, lambda1, lambda2)), digits = 2)

with(mdata, hist(Y, prob = TRUE, main = "Orange=estimate, blue=truth"))
abline(v = 1 / Coef(fit)[c(2, 3)],  lty = 2, col = "orange", lwd = 2)
abline(v = 1 / c(lambda1, lambda2), lty = 2, col = "blue", lwd = 2)

## End(Not run)

###NOT WORKING 

#a function to compare two species logits
logitPlot <- function(dataML_tall, focsp1, focsp2){
  pfocsp1 <- NULL 
  pfocsp2 <- NULL
  #get dataML (see appropriate 'output' ... )
  dataML_tall %>%
    filter(Date_Time =="2020-06-12 05:30:00") %>%
    filter(point == "408") %>%
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
    filter(Date_Time =="2020-06-12 05:30:00") %>%
    filter(point == "408") %>%
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
