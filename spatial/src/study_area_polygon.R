#study area polygon
#study_area_polygon.R
#avian points + buffer
#veg points + buffer for reference
library(sf)
library(mapview)
library(RColorBrewer)
library(here)
library(tidyverse)

#depends upon the total points from "both_points_v" 
#generated in "read_points_output_data.R"
plot <- both_points_v #creates a local object so we don't mess up the original

#encode point type for the "plot" description
plot %>% 
  mutate(type = case_when(veg_point != 0 & avian_point != 0 ~ "both",
                                veg_point == 0  & avian_point != 0 ~ "avian only",
                                veg_point != 0 & avian_point == 0 ~ "veg. only")
  ) -> plot

#bird only
plot %>% 
  dplyr::select(avian_point,
                geometry.vegetation)  %>%
  filter(avian_point != 0) -> df.avian.sf #bird sf data 

#veg only
plot %>% 
  dplyr::select(veg_point,
                geometry.vegetation)  %>%
  filter(veg_point != 0) -> df.vegetation.sf #veg sf data

#using sf::st_convex_hull(x) make a buffered minimum convex polygon for birds
avian_hull <- df.avian.sf  %>%
  summarise(geometry.grp.bird.pt = st_combine(geometry.vegetation))  %>%
  st_convex_hull() %>%
  st_buffer(dist = 400,
            endCapStyle = "ROUND",
            joinStyle = "ROUND")

#... and plants
veg_hull <- df.vegetation.sf  %>%
  summarise(geometry.grp.veg.pt = st_combine(geometry.vegetation))  %>%
  st_convex_hull() %>%
  st_buffer(dist = 400,
            endCapStyle = "ROUND",
            joinStyle = "ROUND")  

#view results
mapView(cbi, alpha.regions = 0.8) +
mapView(cbi_cuts) +
# mapView(cbi4) +
# mapView(cbi4_cuts) +
# mapView(vegetation_points) +
mapView(avian_points, ) +
# mapView(plot, 
#         zcol = "type", 
#         col.regions=brewer.pal(3, "YlGn"), 
#         label = paste0("bird ", df.avian.sf$avian_point,", veg ", df.vegetation.sf$veg_point)) +
mapView(list(avian_hull, veg_hull),
        layer.name = c("bird","vegetation"), 
        col.regions=brewer.pal(3, "BrBG")[c(1,3)])

mapView(list(avian_hull, veg_hull),
        layer.name = c("bird","vegetation"), 
        col.regions=brewer.pal(3, "BrBG")[c(1,3)]) +
mapView(avian_points, ) +
mapView(cbi_cuts) +
mapView(fire_boundary)

exactextractr::exact_extract(cbi_cuts, avian_hull, c("mean","median","min","max","count","sum"))

#stats for study area paragraph
sf::st_intersection(fire_boundary, avian_hull, cbi_cuts) -> test_intersection

exactextractr::exact_extract(cbi_cuts, test_intersection, c("mean","median","min","max","count","sum"))

# find overall min and max extents x ...
min_ext_x <- min(extent(fire_boundary)[1],extent(avian_hull)[1], extent(cbi_cuts)[1])
max_ext_x <- max(extent(fire_boundary)[2],extent(avian_hull)[2], extent(cbi_cuts)[2])
# and y ... 
min_ext_y <- min(extent(fire_boundary)[3],extent(avian_hull)[3], extent(cbi_cuts)[3])
max_ext_y <- max(extent(fire_boundary)[4],extent(avian_hull)[4], extent(cbi_cuts)[4])

#putting it together
overall_extent <- extent(fire_boundary)
overall_extent[1] <- min(min_ext_x, extent(fire_boundary)[1])
overall_extent[2] <- max_ext_x
overall_extent[3] <- min(min_ext_y, extent(fire_boundary)[3])
overall_extent[4] <- max_ext_y

sf::st_crop(fire_boundary, overall_extent)

extent(overall_extent)
extent(fire_boundary)

#save results to output 
#[X] avian_hull
sf::st_write(obj = avian_hull, 
             dsn = here("spatial", "output", "objects","avian_hull.shp"))


sf::st_write(obj = avian_hull, 
             dsn = here("spatial", "output", "objects","avian_hull.csv"), 
             layer_options = "GEOMETRY=AS_WKT",
             delete_dsn = TRUE)

#[X] veg_hull
sf::st_write(obj = veg_hull, 
             dsn = here("spatial", "output", "objects","veg_hull.csv"), 
             layer_options = "GEOMETRY=AS_WKT",
             delete_dsn = TRUE)
