#large tree demo

#data on tree height extracted from raster for 1 ha, one point == canopy_ht_2018_1ha
#an image of canopy height
par(pty="s")
par(mfrow = c(2, 2))
par(mar = c(2,2,1,1))
canopy_ht_2018_1ha[[3]][,1] %>% 
  matrix(.,11,11,byrow=F) %>% 
  raster::image(., xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  title(main = "canopy height")

#large trees defined by being greater or equal to 22m
canopy_ht_2018_1ha[[3]][,1] %>% 
  matrix(.,11,11,byrow=F) >= 22 -> large_tree
  large_tree <- large_tree * 1 # make a number
  raster::image(large_tree, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  title(main = "large tree (canopy height >= 22m")
  
#coverage fraction -- on the edges isn't always 1
canopy_ht_2018_1ha[[3]][,2] %>%   # the second column is the coverage fraction
    matrix(.,11,11,byrow=F) -> coverage_fraction  
    raster::image(coverage_fraction, xaxt= "n", yaxt= "n") %>% 
    axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
    axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
    title(main = "area coverage fraction")

#plot of the values multiplied together  
large_tree * coverage_fraction -> large_trees_coverage_fraction 
  raster::image(large_trees_coverage_fraction, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  title(main = "large tree * coverage fraction")

#the value of the percent of the pixels covered by large trees is 
perc_large_tree <- sum(large_trees_coverage_fraction)/sum(coverage_fraction)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("% large tree =", round(perc_large_tree,3)))

#what is the cover of the large trees?  IS THIS A MEANINGFUL STATISTIC?


canopy_cov_2018_1ha[[3]][,1] %>% 
  matrix(.,11,11,byrow=F) -> canopy_cover
  
  
#coverage fraction though on the edges isn't always 1
canopy_cov_2018_1ha[[3]][,2] %>%   # the second column is the coverage fraction
  matrix(.,11,11,byrow=F)  == coverage_fraction#its the same

large_trees_canopy_cover <- large_tree * canopy_cover * coverage_fraction

#plot of these
par(mfrow = c(1, 3))
par(mar = c(2,2,2,1))
raster::image(large_trees_coverage_fraction, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  title(main = "large tree")

raster::image(canopy_cover, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA)
  title(main = "canopy cover")

raster::image(large_trees_canopy_cover, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(0,10,by=1),xpd=NA)
  title(main = "large_tree * canopy_cover \n* coverage_fraction")
  
#the value of the percent of the pixels covered by large trees is 
cov_large_tree <- sum(large_trees_canopy_cover)/sum(canopy_cover)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("fraction of canopy cover due to large trees =", round(cov_large_tree,3)))

#but the average canopy cover of large trees is simply
mean(large_trees_canopy_cover[large_trees_canopy_cover>0]) -> meanltcov
#[1] 42.34715

text(.5,.3, paste("mean canopy cover of large trees =", round(meanltcov,3)))

