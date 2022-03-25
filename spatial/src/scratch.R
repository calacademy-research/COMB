#random non-cleaned up scratch work ... ignore and don't leave good stuff here

data.rot.melt <- reshape2::melt(data.pca$rotation[,1:5]) # Melt matrix (needed for ggplot), factors 1-5
ggplot(data = data.rot.melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0.0, limit = c(-0.6,0.6),   name="Factor loadings") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()



# Here we’ll show how to calculate the PCA results for variables: coordinates, cos2 and contributions:
#   
#   var.coord = loadings * the component standard deviations
# var.cos2 = var.coord^2
# var.contrib. The contribution of a variable to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)












# Compute the mean direction of a random sample of observations.
x <- circular(runif(50, circular(0), pi))
mean(x)



# Compute summary statistics of a random sample of observations. 
data <- circular(runif(50, 0, pi))
summary(data)
summary(data.frame(data, runif(50, 0, pi)))

x <- rvonmises(n=10, mu=circular(0), kappa=9, control.circular=list(units="degrees"))
par(mfcol=c(2, 2))
plot(x)
y <- conversion.circular(x) # only the unit is changed (to radians) and 
####### the data converted.
plot(y)
z <- conversion.circular(x, units="degrees", zero=pi) # only the zero is changed and 
####### the data converted.
plot(z)
w <- conversion.circular(x, zero=pi/2, rotation="clock") # zero and rotation is 
####### changed and the data converted.
plot(w)


#test
compass<-circular(c(0,90,180,270,360), units = "degrees")
compass_deg<-conversion.circular(compass, units="degrees", zero=pi/2, rotation="clock")
compass_rad<-conversion.circular(compass, units="radians", zero=pi/2, rotation="clock")
plot(compass_deg)
plot(compass_rad)

con_circle<-function(degs) {
  conversion.circular(
    circular(degs, units = "degrees"), 
    units="radians", zero=pi/2, rotation="clock")
  }

# Compute the mean direction of a random sample of observations.
x <- circular(runif(50, circular(0), pi))
mean(x)

plot(x)
points(mean(x), col='red')
# 
# lapply(extract_dem, function(tax){
#   if (!is.null(tax)) tax %>% stack() %>% dplyr::rename(source_taxon_name = "values", source_taxon_rank = "ind")
# }) 

# Compute summary statistics of a random sample of observations. 
data <- circular(runif(50, 0, pi))
summary(data)
summary(data.frame(data, runif(50, 0, pi)))


l1 <- list(list(a = 1L), list(a = NULL, b = 2L), list(b = 3L))
l1 %>% map("a", .default = "???")


x <- rvonmises(n=100, mu=circular(pi), kappa=2)
res25 <- density(x, bw=25, control.circular=list(units="degrees"))
circularp(res25$x)
plot(res25, points.plot=TRUE, xlim=c(-1.6,1))
res50 <- density(x, bw=25, adjust=2)
lines(res50, col=2)
lines(res50, col=3, shrink=0.9) #shrink the plot wrt the function :-)
lines(res50, col=4, offset=0.5) #draw it with a reference circle of 0.5

#to convert degrees on a compass to radians on a compass
par(mfrow=c(1,2))
tstdeg<-1:3.6*100
tstdeg<-circular(tstdeg,units='degrees', zero=pi/2,rotation="clock")
plot(tst,cex=tst/100)
tstrad<-conversion.circular(tstdeg,units='radians') #specificity not required b/c zero and rotation inherited
plot(tstrad, cex=tst/100)


plot(circular(wide_forest_variables_mean_perc_bin_1ha$mean_aspCirc_1ha, units = "radians", zero=pi/2, rotation="clock"))
points(mean.circular(wide_forest_variables_mean_perc_bin_1ha$mean_aspCirc_1ha,units = "radians", zero=pi/2, rotation="clock"), col="red")

par(mfrow=c(2,2))

plot(circular(wide_forest_variables_mean_perc_bin_4ha$mean_aspCirc_4ha, units = "radians", zero=pi/2, rotation="clock"))

plot(circular(wide_forest_variables_mean_perc_bin_4ha$mean_aspCirc_4ha, units = "radians", zero=pi/2, rotation="clock"))
points(mean.circular(wide_forest_variables_mean_perc_bin_1ha$mean_aspCirc_1ha,units = "degrees", zero=pi/2, rotation="clock"), col="red")
