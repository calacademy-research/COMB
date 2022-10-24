# build_pixel_level_data.R
# Pixel by pixel comparisons around particular points 2018 vs 2020
# Depends on canopy_imgStack
# install.packages("naturalsort")
library(naturalsort)

#RAVG addition to wldf_100 (bootstrap style)
#for wldf_100 see read_points_output_data.R line 202
#generate data for Caples RAVG
#see weighted mean function wmf 
exactextractr::exact_extract(RAVG_Caples_imgStack$ca3872412014620191010_20181118_20191118_rdnbr_cbi4, wldf_100, wmf) -> RAVG_1819

RAVG_pct_cutoffs <-c(0, .25, .75, 1)

RAVG_q_cuts <- quantile(RAVG_1819[RAVG_1819 > 0], RAVG_pct_cutoffs)
#note the hack to make it ignore zero

#fix first to equal zero
RAVG_q_cuts[1] <- 0

RAVG_1819 <- cbind(RAVG_1819, cut(RAVG_1819, RAVG_q_cuts)) 
colnames(RAVG_1819)[2] <- c("RAVG_category")

RAVG_1819 %>%
  as_tibble() %>%
  mutate(RAVG_category = ifelse(is.na(RAVG_category), 0, RAVG_category)) -> RAVG_1819
#this finalizes our RAVG categorical variable to be 0,1,2,3 for none, low, moderate, and high severity.
wldf_100 <- cbind(wldf_100, RAVG_1819)

#generate data for CH
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2018, wldf_100) -> test_ch18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2020, wldf_100) -> test_ch20 # %>% head()

#and CC
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2018, wldf_100) -> test_cc18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2020, wldf_100) -> test_cc20 # %>% head()

#add this metadata when binding "pixel_level_db_ch_cc" together

tibble::as_tibble(test_ch18, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100$avian_point, 441)) %>%
  mutate(RAVG_1819 = rep(wldf_100$RAVG_1819,441)) %>%
  mutate(RAVG_category = rep(wldf_100$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_height_18 = value) %>% #colnames()
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_ch18_long

tibble::as_tibble(test_ch20, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100$avian_point, 441)) %>%
  mutate(RAVG_1819 = rep(wldf_100$RAVG_1819,441)) %>%
  mutate(RAVG_category = rep(wldf_100$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_height_20 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_ch20_long

tibble::as_tibble(test_cc18, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100$avian_point, 441)) %>%
  mutate(RAVG_1819 = rep(wldf_100$RAVG_1819,441)) %>%
  mutate(RAVG_category = rep(wldf_100$RAVG_category,441)) %>%
  mutate(point_idx = name) %>% 
  mutate(canopy_cover_18 = value) %>% 
  select(c(9,4:8,10)) %>% 
  arrange(point_idx) -> test_cc18_long

tibble::as_tibble(test_cc20, .rows = 441, .name_repair = "universal") %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(cols = -rowname) %>% 
  mutate(pixel_name = rowname) %>% 
  mutate(veg_point = rep(wldf_100$veg_point, 441)) %>% 
  mutate(avian_point = rep(wldf_100$avian_point, 441)) %>%
  mutate(RAVG_1819 = rep(wldf_100$RAVG_1819,441)) %>%
  mutate(RAVG_category = rep(wldf_100$RAVG_category,441)) %>%
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
                                     RAVG_1819, 
                                     .fun = min, 
                                     .desc = FALSE)) -> pixel_level_db_ch_cc

jitter <- position_jitter(width = 0.15, height = 0.15)
  
    pixel_level_db_ch_cc %>%
    # filter(avian_point != 0) %>%
    filter(avian_point %in% c(841, 454, 490,576, 1057, 1072)) %>%
    ggplot(aes(x = canopy_height_18, y = canopy_height_20)) +
    # geom_point(pch = ".") +
    # stat_summary(fun.data=mean_cl_normal) + 
    # geom_point(aes(size = coverage_fraction, alpha = 1)) +
    geom_point(aes(size = coverage_fraction), alpha = 0.42, position = jitter, shape = 21, fill = "tan", stroke = .2) +
    # geom_point(position = jitter) +
    scale_size(range = c(0,3)) +
      geom_smooth(method="lm", 
                formula= y~x,
                mapping = aes(weight = coverage_fraction),
                color = "red",
                lwd = 0.5,
                show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = "blue") +
    coord_fixed() +
      xlim(0, 40) +
      ylim(0, 40) +
    theme(aspect.ratio=1) +
    theme(legend.position="none") +
    facet_wrap(~avian_point) # +
    # theme(strip.text = element_blank(), panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))

  
    # scale_color_manual(values = c("high" = "black",
    #                               "moderate" = "red",
    #                               "low" = "orange",
    #                               "unburned" = "green")) +
    # ggtitle(c("Relationship between canopy variables vs \nburn severity (RAVG) before and after fire")) +
    # xlab("Change in height (2020-2018)") +
    # ylab("Change in cover (2020-2018)") +
    # facet_wrap(~RAVG) 
  