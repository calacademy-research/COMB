#large tree demo

#data on tree height extracted from raster for 1 ha, one point == canopy_ht_2018_1ha
#an image of canopy height
par(pty="s")
par(mfrow = c(2, 2))
par(mar = c(1,2,1,1))
canopy_ht_2018_1ha[[3]][,1] -> canopy_ht 
  matrix(canopy_ht,11,11,byrow=F) %>% 
  raster::image(., xaxt= "n", yaxt= "n",xlab = "", ylab="") %>% 
    axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
    axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
    title(main = "canopy height")

#the mean canopy height is 
mean_canopy_height <- sum(canopy_ht)/sum(coverage_fraction)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("mean canopy height =", round(mean_canopy_height,3)))


#large trees defined by being greater or equal to 22m
canopy_ht_2018_1ha[[3]][,1] %>% 
  matrix(.,11,11,byrow=F) >= 22 -> large_tree
  large_tree <- large_tree * 1 # make a number
  raster::image(large_tree, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  title(main = "large tree (canopy height >= 22m")
#the mean canopy height of large trees is 
  mean_canopy_height_lt <- sum(canopy_ht*large_tree)/sum(coverage_fraction*large_tree)
  #note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
  #putting the value in there
  text(.5,.5, paste("mean canopy height of large trees =", round(mean_canopy_height_lt,3)))
  
#coverage fraction -- on the edges isn't always 1
canopy_ht_2018_1ha[[3]][,2] %>%   # the second column is the coverage fraction
    matrix(.,11,11,byrow=F) -> coverage_fraction  
    raster::image(coverage_fraction, xaxt= "n", yaxt= "n") %>% 
      axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
      axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
      title(main = "area coverage fraction")
#the mean canopy height of large trees is 
  sum_coverage_fraction <- sum(coverage_fraction)
    #note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
    #putting the value in there
    text(.5,.5, paste("sum coverage fraction =", round(sum_coverage_fraction,3)))
    
#plot of the values multiplied together large tree coverage fraction
large_tree * coverage_fraction -> large_trees_coverage_fraction 
  raster::image(large_trees_coverage_fraction, xaxt= "n", yaxt= "n") %>% 
    axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
    axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
    title(main = "large tree * coverage fraction")

#the value of the percent of the pixels covered by large trees is 
perc_large_tree <- sum(large_trees_coverage_fraction)/sum(coverage_fraction)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("% large tree =", round(perc_large_tree,3)))


#question #2
#what is the cover of the large trees?  IS THIS A MEANINGFUL STATISTIC?
canopy_cov_2018_1ha[[3]][,1] %>% 
  matrix(.,11,11,byrow=F) -> canopy_cover
  
  
#coverage fraction though on the edges isn't always 1
canopy_cov_2018_1ha[[3]][,2] %>%   # the second column is the coverage fraction
  matrix(.,11,11,byrow=F)  == coverage_fraction#its the same

large_trees_canopy_cover <- large_tree * canopy_cover * coverage_fraction

#plot of these
#large_trees_coverage_fraction
par(mfrow = c(3, 1))
par(mar = c(2,2,2,1))
raster::image(large_trees_coverage_fraction, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  title(main = "large tree")
#the value of the percent of the pixels covered by large trees is 
perc_large_tree <- sum(large_trees_coverage_fraction)/sum(coverage_fraction)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("% large tree =", round(perc_large_tree,3)))

#canopy cover
raster::image(canopy_cover, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  title(main = "canopy cover")
#the value of the percent of the pixels covered by large trees is 
  canopy_cov_mean <- sum(canopy_cover)/sum(coverage_fraction)
  #note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
  #putting the value in there
  text(.5,.5, paste("mean canopy cover overall =", round(canopy_cov_mean,3)))
  
#large_trees_canopy_cover
raster::image(large_trees_canopy_cover, xaxt= "n", yaxt= "n") %>% 
  axis(side=1,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  axis(side=2,at=seq(0,1,by=.1),labels=seq(1,11,by=1),xpd=NA) %>%
  title(main = "large_tree * canopy_cover \n* coverage_fraction")
  
#the value of the percent of the pixels covered by large trees is 
cov_large_tree <- sum(large_trees_canopy_cover)/sum(canopy_cover)
#note the latter is 100 which is exactly 100 m^2 (with 10, 10 meter pixels on a side)!
#putting the value in there
text(.5,.5, paste("fraction of canopy cover\ndue to large trees =", round(cov_large_tree,3)))

#but the average canopy cover of large trees is simply
mean(large_trees_canopy_cover[large_trees_canopy_cover>0]) -> meanltcov
#[1] 42.34715

text(.5,.3, paste("mean canopy cover of large trees =", round(meanltcov,3)))

