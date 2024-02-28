library(ggmap)
library(dplyr)
library(raster) #masks dplyr::select()

load(file = "caples_map2.RData")
tab <- read_csv("models/output/bbwo20.csv")
load("BBWO2020_allmodels.RData")

pc.psi <- as.data.frame(bbwo20pc$summary[9:89,])
#pc.psi$z <- round(z$mean[1:81], 2)

tab$pc.psi <- pc.psi$mean
tab$z <- round(zmat[1:81,1], 2)

tif <- stack("ca3872412014620191010_20181118_20191118_rdnbr_cbi4.tif")
tif_df <-
  as.data.frame(tif, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit()
head(tif_df)

ggmap(caples_sat) + 
  geom_point(data=tab,
             aes(x = long, y = lat, fill = z),
             size = 5, 
             pch = 21,
             colour = "black") + 
  #guides(fill = "legend")  + 
  geom_point(data= subset(tab, both==0), 
             aes(x = long, y = lat),
             size = 5,
             pch = 4,
             color = "black") + 
  labs(title = "Predicted Occupancy (z): Full Model",
       x = "",
       y = "", 
       fill = "Mean Posterior\nEstimate of z") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_continuous(type="viridis") +
  theme(legend.position = "none", legend.key = element_blank(), 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

#ggmap(caples_sat) + 
#  geom_raster(data=tif_df, aes(x = x, y = y, fill = ca3872412014620191010_20181118_20191118_rdnbr_cbi4), alpha=0.2) +
trans.points <- ggplot() +
    geom_point(data=tab,
             aes(x = long, y = lat, fill = z),
             size = 5, 
             pch = 21,
              colour = "black") + 
  #guides(fill = "legend")  + 
  geom_point(data= subset(tab, both==0), 
             aes(x = long, y = lat),
             size = 5,
             pch = 4,
             color = "black") + 
  labs(title = "Predicted Occupancy: Full Model",
       x = "",
       y = "", 
       fill = "Mean Posterior\nEstimate of z") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_continuous() +
  theme(legend.position = "none", legend.key = element_blank(), 
            panel.background = element_rect(fill='transparent'), #transparent panel bg
            plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank(), #remove minor gridlines
            legend.background = element_rect(fill='transparent'), #transparent legend bg
            legend.box.background = element_rect(fill='transparent') #transparent legend panel
          )

ggsave(plot=trans.points, filename="models/figures/transparent_pts_post_z.png", bg="transparent",
       width = 6, height = 5, units = "in", dpi = 300)


# point count only estimates of psi
ggmap(caples_sat) + 
  #  geom_raster(data=tif_df, aes(x = x, y = y, fill = ca3872412014620191010_20181118_20191118_rdnbr_cbi4), alpha=0.2) +
  geom_point(data = tab, 
             aes(x = long, y = lat, fill = pc.psi),
             size = 3, 
             pch = 21,
             colour = "white") + 
 # guides(fill = "legend")  + 
  geom_point(data= subset(tab, both ==0), 
             aes(x = long, y = lat),
             size = 3,
             pch = 4,
             color = "grey70") + 
  labs(title = "Predicted Occupancy: PC Only",
       x = "",
       y = "", 
       fill = "Mean Posterior\nProbability of \nOccupancy") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_continuous(type = "viridis") +
  theme(legend.position = "bottom", legend.key = element_blank())

# naive occupancy
df <- tab %>% mutate(occ_cat = case_when(
  naive.pc.pa == 1 & naive.aru.pa == 1 ~ "ARU & PC Detected",
  naive.pc.pa == 1 & naive.aru.pa == 0 ~ "PC Only Detected",
  naive.pc.pa == 0 & naive.aru.pa == 1 ~ "ARU Only Detected",
  naive.pc.pa == 0 & naive.aru.pa == 0 ~ "Undetected via both"
))


ggmap(caples_sat) + 
  geom_point(data = df, 
             aes(x = long, y = lat, fill = occ_cat),
             size = 3,
             pch = 21,
             colour = "white") + 
  geom_point(data= subset(df, occ_cat == "Unoccupied via both"), 
             aes(x = long, y = lat),
             size = 3,
             pch = 4,
             color = "grey70")  + 
  labs(title = "Naive Occupancy",
       x = "",
       y = "", 
       fill = "Naive Occupancy") + 
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values = c('#FDE725FF', '#55C667FF','#33638DFF','#440154FF')) +
  theme(legend.position = "bottom", legend.key = element_blank())


