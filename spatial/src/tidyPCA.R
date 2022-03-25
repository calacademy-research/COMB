#pca tidy style 
library(tidyverse)
#devtools::install_github("tidymodels/broom")
#library("broom", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
library(broom)
library(cowplot)

#getting our data together (depends on Caples_Spatial_EDA.Rmd) to make "canopy_nbr_dem_cir"
exactextractr::exact_extract(canopy_fuel_nbr_dem_cir, wldf_100) %>%
  do.call(rbind, .) -> pointdata

# pca_fit <- pointdata %>% 
#   select(-aspect, -coverage_fraction) %>%
#   dplyr::select(where(is.numeric)) %>% # retain only numeric columns
#   prcomp(scale = TRUE) # do PCA on scaled data

pointdata %>% 
  dplyr::select(-aspect, -coverage_fraction, -Caples_dNBR_Nov18_Nov19) %>% #head()
  mutate(aspect_radians=as.numeric(cut(aspect_radians,4))) %>% 
  dplyr::select(where(is.numeric)) %>% # retain only numeric columns
  scale() -> scaledpointdata # scale data

pca_fit <- prcomp(scaledpointdata) # do PCA

#PVE and scree plots
pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

#plot heatmap of factor loadings
varlabs<-c("elevation",
           "Height 2020", "Height 2019", "Height 2018", 
           "Cover 2020", "Cover 2019", "Cover 2018", 
           "NBR Nov20", "NBR Nov19", "NBR Nov18", 
           "aspect binary", "aspect radians")

#get rotation matrix filter out number of components desired
pca_fit %>%
  tidy(matrix = "rotation") %>%
  filter(PC <= 6) -> data.rot
#make the plot
ggplot(data = data.rot, aes(PC, column, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0.0, limit = range(data.rot$value),   name="Factor loadings") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 42, vjust = 1, 
                                   size = 12, hjust = 1)) +
  scale_x_continuous("PC", labels = as.character(1:6), breaks = 1:6) +
  scale_y_discrete("Variables") +# ,labels = varlabs)
  coord_fixed()

#table of summary percent contribution of each variable to each component
# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- pca_fit %>%
  tidy(matrix = "loadings") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  select(-column)
sdev <- pca_fit$sdev
var.coord <- t(apply(as.matrix(loadings), 1, var_coord_func, comp.sdev=sdev)) 
var.coord[, 1:6]

# Compute Cos2
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])

# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
pca_fit %>%
  tidy(matrix = "loadings") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  select(column) -> vars
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
var.contrib<-cbind(vars,var.contrib)
view(var.contrib[, 1:7])

pca_fit %>%
  broom::augment(pointdata) -> aug_pca_fit 

#manually build all combinations of the 'pfit' plots
#should NOT do by hand, but did for 1-2 through 5-6 component combinations
aug_pca_fit %>% #head() # add original dataset back in
  ggplot(aes(x = .fittedPC5, y = .fittedPC6, color = Caples_dNBR_Nov18_Nov19, shape = cut(aspect_radians,4))) + #, color = outcome)) + fix for burn/no burn
  geom_point(size = 1.5) +
  #scale_color_manual(values = c(.fittedPC1 = "#D55E00", .fittedPC2 = "#0072B2")) +
  #scale_color_manual(values = c(Caples_dNBR_Nov18_Nov19 = "#D55E00", aspect_radians = "#0072B2")) +
  theme_half_open(12) + background_grid() -> pfit56

#combined plots of fitted values and rotation matrices 
# for combinations up to 5 & 6.
# define arrow style for plotting
arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)

# plot rotation matrix
# and the rotation matrix plots.
pca_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  ggplot(aes(PC5, PC6)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  geom_text(
    aes(label = column),
    hjust = 1, nudge_x = -0.02, 
    color = "#904C2F"
  ) +
  xlim(-1.25, .5) + ylim(-.5, 1) +
  coord_fixed() + # fix aspect ratio to 1:1
  theme_minimal_grid(12) -> pr56

#build up plots for the whole story
require(gridExtra)
plot1 <- pfit13
plot2 <- pr13
grid.arrange(plot1, plot2, ncol=2)

pca_fit %>%
  tidy(matrix = "eigenvalues")



pca_fit %>%
  unnest(pc_frac) %>%
  knitr::kable()
