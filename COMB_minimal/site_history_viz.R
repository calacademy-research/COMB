# requires raw data for all woodpeckers
# under construction! need loop for species


aru_dets <- read_csv("models/input/file_logit_agg.csv") %>% filter(species=="BBWO", year(Date_Time) == 2020)

# compare with point count detections
pc_dets <- read_csv("../COMB/point_counts/data_ingest/output/PC_delinted.csv") %>% filter(year==2020 & birdCode_fk=="BBWO")

pc_dets$Date <- as.Date(pc_dets$DateTime)
pc_dets$point <- pc_dets$point_ID_fk

head(pc_dets)
pc_dets$fill <- ifelse(pc_dets$abun > 0, "filled", "empty")
sumabun <- pc_dets %>% group_by(point) %>% summarise(sumabun=sum(abun))
pc_dets <- left_join(pc_dets, sumabun)



sitehist_plot <- ggplot() +
  geom_point(data=aru_dets[!is.na(aru_dets$point),],
             aes(x=Date_Time, y=max_logit, color=-max_logit)) +
  geom_point(data=pc_dets,
             aes(x=DateTime, y=abun, fill=fill), shape=21, size=2) +
  geom_hline(yintercept=0) +
  scale_fill_manual(values = c("filled"="black", "empty"="white")) +
  scale_color_viridis_c(option = "magma")+
  facet_wrap(~point) +
  labs(y="logit", x="date", title=aru_dets$species[[1]]) +
  theme_bw() +
  xlim(c(min(aru_dets$Date_Time), max(aru_dets$Date_Time))) +
  ylim(c(-3,10)) +
  theme(legend.position = "none",
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(plot = sitehist_plot, filename=paste("COMB_minimal/results/figures/site_histories/site_history", aru_dets$species[[1]], "2020", Sys.Date(), "plot.png", sep="_"),width = 10, height = 6, units = "in", dpi = 300)


