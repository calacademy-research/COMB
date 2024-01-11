library(ggmap)
library(dplyr)

load(file = "caples_map2.RData")

df <- read.csv('bbwo20 - bbwo20.csv', 
               header = TRUE)

ggmap(caples_sat) + 
  geom_point(data = df, 
             aes(x = long, y = lat, fill = psi),
             size = 5, 
             pch = 21,
              colour = "white") + 
  guides(fill = "legend")  + 
  geom_point(data= subset(df, occ_cat == "Unoccupied via both"), 
             aes(x = long, y = lat),
             size = 5,
             pch = 4,
             color = "grey70") + 
  labs(title = "BBWO 2020 Predicted Occupancy",
       x = "",
       y = "", 
       fill = "Mean Posterior\nProbability of \nOccupancy") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_continuous(type = "viridis")

df <- df %>% mutate(occ_cat = case_when(
  naive.pc.pa == 1 & naive.aru.pa == 1 ~ "ARU & PC Occupied",
  naive.pc.pa == 1 & naive.aru.pa == 0 ~ "PC Only Occupied",
  naive.pc.pa == 0 & naive.aru.pa == 1 ~ "ARU Only Occupied",
  naive.pc.pa == 0 & naive.aru.pa == 0 ~ "Unoccupied via both"
))

ggmap(caples_sat) + 
  geom_point(data = df, 
             aes(x = long, y = lat, fill = factor(occ_cat)),
             size = 5,
             pch = 21,
             colour = "white") + 
  geom_point(data= subset(df, occ_cat == "Unoccupied via both"), 
             aes(x = long, y = lat),
             size = 5,
             pch = 4,
             color = "grey70")  + 
  labs(title = "BBWO 2020 Naive Occupancy",
       x = "",
       y = "", 
       fill = "Naive Occupancy") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values = c('#55C667FF', '#FDE725FF','#33638DFF','#440154FF'))
