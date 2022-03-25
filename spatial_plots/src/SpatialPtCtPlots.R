## Read in the data

year_pt_counts <- read_csv("~/....csv")

## Generating a basic point plot of the sites on top of a Google map
library(RgoogleMaps)
comb_map <- GetMap.bbox(lonR = range(year_pt_counts$longitude, na.rm = TRUE),
                        latR = range(year_pt_counts$latitude, na.rm = TRUE), 
                        size = c(640, 640),
                        maptype = "satellite")

PlotOnStaticMap(comb_map)

convert_points <- LatLon2XY.centered(comb_map,
                                     year_pt_counts$latitude,
                                     year_pt_counts$longitude)
points(convert_points$newX, convert_points$newY, col = 'red', pch=19)

## Plotting the point counts as the character size for bird = "hamfly"
library(spBayes)
library(classInt)
library(RColorBrewer)

hamfly_df <- subset(subset(year_pt_counts, (bird == "hamfly") & (year == 2018)), !is.na(latitude))
hamfly_ptct <- hamfly_df$detection_count
hamfly_maxlgt <- hamfly_df$max_logit

coords <- as.matrix(hamfly_df[, c("longitude", "latitude")])
plot(coords, pch = 1, cex = sqrt(hamfly_ptct)/10, col = "darkgreen")

# Generating interpolation plots
library(MBA)
library(fields)
x.res <- 100
y.res <- 100
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
## Point count surface plot
surf <- mba.surf(cbind(coords, hamfly_ptct), no.X = x.res, no.Y = y.res, 
                 h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Lon", ylab = "Lat", col = col.br(25), main = "Point Count - Hammonds Flycatcher")
contour(surf, add=TRUE)

## Max logit surface plot
surf <- mba.surf(cbind(coords, hamfly_maxlgt), no.X = x.res, no.Y = y.res, 
                 h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Lon", ylab = "Lat", col = col.br(25), main = "Max Logit - Hammonds Flycatcher")
contour(surf, add=TRUE)