focsp1 <- "NOGO"
dataML_top1 %>%
  ungroup() %>% # as_tibble() %>%
  # filter(!is.na(point)) %>% dim()
  slice_sample(.data = ., n = 10^7) %>% 
  #  filter(Date_Time =="2020-06-11 09:30:00") %>%
  # filter(point == "444") %>%
  arrange(Start_Time) %>% 
  mutate(focal_species = if_else(species == focsp1, "yes", "no")) %>%
  ggplot(., aes(x=logit, color = focal_species)) +
  geom_density(adjust = 1) +
  xlim(-3, 4) +
  labs(title=paste("Logits for a single species",
                   focsp1, "versus all others")) +
  # Add mean line
  geom_vline(data = . %>%
               group_by(focal_species) %>%
               summarize(focal_species_mean = mean(logit)),
             mapping = aes(xintercept = focal_species_mean),
             color = c(2,3), linetype="dashed", size=1) -> pfocsp1

HAWO_t1<- pfocsp1
BUSH_t1<- pfocsp1
