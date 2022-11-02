#study area polygon
#study_area_polygon.R
#avian points + buffer
#veg points + buffer for reference
library(sf)
library(mapview)
library(RColorBrewer)
library(here)

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
mapView(plot, 
        zcol = "type", 
        col.regions=brewer.pal(3, "YlGn"), 
        label = paste0("bird ", df.sf$avian_point,", veg ", df.sf$veg_point)) +
mapView(list(avian_hull, veg_hull),
        layer.name = c("bird","vegetation"), 
        col.regions=brewer.pal(3, "BrBG")[c(1,3)])

#save results to output 
#[X] avian_hull
sf::st_write(obj = avian_hull, 
             dsn = here("spatial", "output", "objects","avian_hull.csv"), 
             layer_options = "GEOMETRY=AS_WKT",
             delete_dsn = TRUE)

#[X] veg_hull
sf::st_write(obj = veg_hull, 
             dsn = here("spatial", "output", "objects","veg_hull.csv"), 
             layer_options = "GEOMETRY=AS_WKT",
             delete_dsn = TRUE)
