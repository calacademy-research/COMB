# build_pixel_level_data.R
# Pixel by pixel comparisons around particular points 2018 vs 2020
# Depends on canopy_imgStack
# install.packages("naturalsort")
library(naturalsort)

#RAVG addition to wldf_100 (bootstrap style)
#for wldf_100 see read_points_output_data.R line 202
#generate data for Caples RAVG
#see weighted mean function wmf 
exactextractr::exact_extract(RAVG_Caples_imgStack$ca3872412014620191010_20181118_20191118_rdnbr_cbi4, wldf_100, wmf) -> RAVG_1819_cbi4

#test new ravg cutoffs (see AW email)

exactextractr::exact_extract(RAVG_Caples_imgStack$ca3872412014620191010_20181118_20191118_rdnbr_cbi, wldf_100, wmf) -> RAVG_1819_cbi


#RAVG_pct_cutoffs <-c(0, .25, .75, 1)

#RAVG_q_cuts <- quantile(RAVG_1819[RAVG_1819 > 0], RAVG_pct_cutoffs)
#note the hack to make it ignore zero

#fix first to equal zero
#RAVG_q_cuts[1] <- 0

#RAVG_f_cuts (f = fixed cuts)
RAVG_f_cuts_cbi4 <- c(-1, 0, 1, 3, 4)
RAVG_f_cuts_cbi4_label = c("unburned 0", "low (0,1]", "moderate (1,3]", "high (3,4]")

#from Angela's email 2022-10-27
#and https://burnseverity.cr.usgs.gov/ravg/background-products-applications
#Outside fire perimeter (0)
#Unchanged (0 <= CBI < 0.1) 
#Low severity (0.1 <= CBI < 1.25) 
#Moderate severity (1.25 <= CBI < 2.25) 
#High severity (2.25 <= CBI < 3) 

#0 < 0.1 < 1.25 < 2.25 < 3
RAVG_f_cuts_cbi_to_cbi4 <- c(-1, 0, .1, 1.25, 2.25, 3)
RAVG_f_cuts_cbi_to_cbi4_label = c("outside burn perimeter 0", "unburned 0 - 0.1", "low (0.1,1.25]", "moderate (1.25,2.25]", "high (2.25,3]")


RAVG_1819_cbi4 <- cbind(as_tibble(RAVG_1819_cbi4), cut(RAVG_1819_cbi4, breaks=RAVG_f_cuts_cbi4, labels = RAVG_f_cuts_cbi4_label)) 
colnames(RAVG_1819_cbi4)[1] <- c("RAVG_1819_cbi4")
colnames(RAVG_1819_cbi4)[2] <- c("RAVG_cbi4_category")

RAVG_1819_cbi <- cbind(as_tibble(RAVG_1819_cbi), cut(RAVG_1819_cbi, breaks=RAVG_f_cuts_cbi_to_cbi4, labels = RAVG_f_cuts_cbi_to_cbi4_label)) 
colnames(RAVG_1819_cbi)[3] <- c("RAVG_1819_cbi")
colnames(RAVG_1819_cbi)[3] <- c("RAVG_cbi_to_cbi4_category")

RAVG_1819 <- cbind(RAVG_1819_cbi4, RAVG_1819_cbi)

RAVG_1819 %>% View()
  # as_tibble() %>%
  # mutate(RAVG_category = ifelse(is.na(RAVG_category), 0, RAVG_category)) -> RAVG_1819
#this finalizes our RAVG categorical variable to be 0,1,2,3 for none, low, moderate, and high severity.
wldf_100_wRAVG <- cbind(wldf_100, RAVG_1819)

#generate data for CH
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2018, wldf_100_wRAVG) -> test_ch18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2020, wldf_100_wRAVG) -> test_ch20 # %>% head()

#and CC
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2018, wldf_100_wRAVG) -> test_cc18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2020, wldf_100_wRAVG) -> test_cc20 # %>% head()

#add this metadata when binding "pixel_level_db_ch_cc" together

tibble::as_tibble(test_ch18, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>%  
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100_wRAVG$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100_wRAVG$avian_point, 441)) %>%
  mutate(RAVG_1819_cbi4 = rep(wldf_100_wRAVG$RAVG_1819_cbi4,441)) %>%
  mutate(RAVG_category = rep(wldf_100_wRAVG$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_height_18 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_ch18_long

tibble::as_tibble(test_ch20, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100_wRAVG$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100_wRAVG$avian_point, 441)) %>%
  mutate(RAVG_1819_cbi4 = rep(wldf_100_wRAVG$RAVG_1819_cbi4,441)) %>%
  mutate(RAVG_category = rep(wldf_100_wRAVG$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_height_20 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_ch20_long

tibble::as_tibble(test_cc18, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100_wRAVG$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100_wRAVG$avian_point, 441)) %>%
  mutate(RAVG_1819_cbi4 = rep(wldf_100_wRAVG$RAVG_1819_cbi4,441)) %>%
  mutate(RAVG_category = rep(wldf_100_wRAVG$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_cover_18 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_cc18_long

tibble::as_tibble(test_cc20, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100_wRAVG$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100_wRAVG$avian_point, 441)) %>%
  mutate(RAVG_1819_cbi4 = rep(wldf_100_wRAVG$RAVG_1819_cbi4,441)) %>%
  mutate(RAVG_category = rep(wldf_100_wRAVG$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_cover_20 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_cc20_long

bind_cols(test_ch18_long,
          test_ch20_long$canopy_height_20,
          test_cc18_long$canopy_cover_18,
          test_cc20_long$canopy_cover_20) %>% 
  unnest(cols = c(canopy_height_18)) %>% 
  select(c(1,3,4,2,5,6,8,7,9,11,13)) -> pixel_level_db_ch_cc
  
  colnames(pixel_level_db_ch_cc)[8:11] <- c("canopy_height_18",
                                      "canopy_height_20",
                                      "canopy_cover_18",
                                      "canopy_cover_20")
  
#fix the levels of the points
#e.g. change the levels(pixel_level_db_ch_cc$avian_point) 
    
pixel_level_db_ch_cc %>%
    mutate(avian_point = as_factor(avian_point)) %>%
    mutate(avian_point = fct_reorder(avian_point, 
                                     RAVG_1819_cbi4, 
                                     .fun = min, 
                                     .desc = FALSE)) -> pixel_level_db_ch_cc

#investigate RAVG differences
par(mfrow=c(1,1))

cbi_to_cbi4_category <- RAVG_1819_cbi4$RAVG_cbi_to_cbi4_category
cbi4_category <- RAVG_1819_cbi4$RAVG_cbi4_category

table(cbi_to_cbi4_category,cbi4_category)

cexval <- .875
plot(y = jitter(RAVG_1819$RAVG_1819_cbi4), 
     x = jitter(RAVG_1819$RAVG_1819_cbi), 
     xlab = "RAVG CBI (continuous)", 
     ylab = "RAVG 'CBI4' (also continuous)")
abline(v=c(0.1,1.25,2.25), lty = 3)
text(x = 0, y = .5, c("unchanged"), srt = 90, cex = cexval, col = "green")
text(x = .1 + (1.25-.1)/2, y = .5, c("low severity"), srt = 90, cex = cexval, col = "brown")
text(x = 1.25 + (2.25-1.25)/2, y = 2, c("moderate severity"), srt = 90, cex = cexval, col = "orange")
text(x = 2.25 + ((3-2.25)/2), y = 3.5, c("high severity"), srt = 90, cex = cexval, col = "red")
abline(h=c(0, 1, 3), lty=2)


df <- cbind(datax,datay)

ggplot(data = df, 
       aes(x=datax, fill=datay)) + 
  geom_bar(position = "fill") + ylab("proportion") #+
  
  
  ggplot(data = interviews_plotting, aes(x = respondent_wall_type, fill = village)) +
  geom_bar(position = "fill") + ylab("proportion") +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white")
  
  
  

plot(x = datax, 
     y = datay, 
     xaxt = "n", 
     xlab = "")

# , 
#      ylab="cbi4 cut to 4 categories", 
#      xlab="cbi cut to cbi4",
#      xaxt = "n", ann=FALSE,
#      yaxt = "n", ann=FALSE,
#      axis(side = 4, labels = FALSE))

# Draw x-axis without labels.


## Draw the x-axis labels.
text(x = cumsum(table(datax)/sum(table(datax))),
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.1,
     ## Use names from the data list.
     labels = levels(RAVG_1819_cbi4$RAVG_cbi_to_cbi4_category),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 1.1,
     ## Increase label size.
     cex = .9)

plot(RAVG_1819_cbi4$RAVG_1819_cbi4, RAVG_1819$RAVG_1819_cbi)

#make a plot of the final story (pick your RAVG)
jitter <- position_jitter(width = 0.15, height = 0.15)

#make mixed effects model 
library(lme4)
model <- lmer(canopy_height_20 ~ canopy_height_18 + RAVG_category + (1|avian_point),
              weight = coverage_fraction, 
              data=pixel_level_db_ch_cc)

pixel_level_db_ch_cc$fit <- predict(model)   #Add model fits to dataframe

pt_colors <- c("dark green","orange","red","black")
  
# plot it
pixel_level_db_ch_cc %>% 
    filter(avian_point != 0) %>%
    # filter(avian_point %in% c(841, 454, 490,576, 1057, 1072)) %>%
    ggplot(aes(x = canopy_cover_18, y = canopy_cover_20)) +
    # geom_point(pch = ".") +
    # stat_summary(fun.data=mean_cl_normal) + 
    geom_point(aes(size = coverage_fraction, col=RAVG_category), 
               alpha = 0.42, 
               position = jitter, 
               shape = 1, 
               # fill = "tan", 
               stroke = .5) +
    scale_color_manual(values=pt_colors) +
    scale_size(range = c(0,1)) +
      geom_smooth(method="lm",
                formula= y~x,
                mapping = aes(weight = coverage_fraction),
                color = "purple",
                lwd = 0.5,
                show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
    geom_vline(xintercept = 22, lty = 3, col = "black") +
    geom_hline(yintercept = 22, lty = 3, col = "black") +
    coord_fixed(ratio = 1, xlim = c(0,100), ylim = c(0,100)) +
      # xlim(0, 40) +
      # ylim(0, 40) +
    # theme(aspect.ratio=1) +
    theme(legend.position="none") +
    facet_wrap(~avian_point) + 
    # geom_vline(aes(xintercept = mean, group = avian_point), colour = 'red') +
    #theme(strip.text = element_text(size = rel(3.0), vjust = -4.0), 
    #     panel.spacing.y = unit(-2, "lines")) # +
     theme(strip.text = element_blank(), 
         panel.spacing.y = unit(0.1, "lines"),
         panel.spacing.x = unit(0.1, "lines")) +
         theme(axis.text.x=element_text(size=7, angle = 90)) +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
         theme(axis.text.y=element_text(size=7)) +
  theme(axis.line = element_line(color='black'),
        # plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

## modeling along the lines of this analysis
### Build full model
### Project 