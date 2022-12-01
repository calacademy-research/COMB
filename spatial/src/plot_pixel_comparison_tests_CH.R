# Pixel by pixel comparisons around particular points 2018 vs 2020
# Depends on canopy_imgStack
install.packages("naturalsort")
library(naturalsort)
#generate data for CH
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2018, wldf_100) -> test_ch18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2020, wldf_100) -> test_ch20 # %>% head()

#and CC
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2018, wldf_100) -> test_cc18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2020, wldf_100) -> test_cc20 # %>% head()

# plot(canopy_imgStack$CaplesCanopyHeight2018, canopy_imgStack$CaplesCanopyHeight2020)
# 
# View(wldf_100)
for(i in naturalsort::naturalsort(wldf_100$avian_point[wldf_100$avian_point != 0])){
# plot 450 veg 106
avian_point_number <- i
rowidno <- as.numeric(dimnames(wldf_100)[[1]][wldf_100$avian_point==avian_point_number])

par(mfrow=c(1,1))
par(pty="s")
par(mar=c(4,4,3,2))
plot(jitter(test_ch18[[rowidno]][,1]), jitter(test_ch20[[rowidno]][,1]), xlab = "Canopy Height in 2018", ylab = "Canopy Height in 2020",pch = ".", cex = 3)
abline(0,1)
abline(h = 22, v = 22)
abline(lm(test_ch20[[rowidno]][,1]~test_ch18[[rowidno]][,1]), col="red")
title(main=paste0("Pixel by pixel comparison \n for point ", avian_point_number))  #, " 2018 vs 2020"

#abline(lm(test_ch20[[rowidno]][,1] ~ 1 + offset(1* test_ch18[[rowidno]][,1])), col="red")

test_ch18[[rowidno]][,1] >= 22 -> ch2018g22
sum(ch2018g22 * 1)/length(test_ch18[[rowidno]][,1]) -> pctCHg22_2018
text(x = max(test_ch20[[rowidno]][,1])*.7, y = max(test_ch18[[rowidno]][,1])*.2,labels = paste(as.character(round(pctCHg22_2018,2)*100), "% \nlarge tree"))

test_ch20[[rowidno]][,1] >= 22 -> ch2020g22
sum(ch2020g22 * 1)/length(test_ch20[[rowidno]][,1])  -> pctCHg22_2020
text(x = max(test_ch18[[rowidno]][,1])*.2, y = max(test_ch20[[rowidno]][,1])*.9,labels = paste(as.character(round(pctCHg22_2020,2)*100), " % \nlarge tree"))

par(mfrow=c(1,1))
par(pty="s")
par(mar=c(4,4,3,2))
plot(jitter(test_cc18[[rowidno]][,1]), jitter(test_cc20[[rowidno]][,1]), xlab = "Canopy Cover in 2018", ylab = "Canopy Cover in 2020", pch = ".", cex = 3)
abline(0,1)
# abline(h = 22, v = 22)
abline(lm(test_cc20[[rowidno]][,1]~test_cc18[[rowidno]][,1]), col="red")
title(main=paste0("Pixel by pixel comparison \n for point ", avian_point_number))  #, " 2018 vs 2020"

}
