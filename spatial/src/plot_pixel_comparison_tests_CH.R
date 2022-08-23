exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2018, wldf_100) -> test_ch18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2020, wldf_100) -> test_ch20 # %>% head()

plot(canopy_imgStack$CaplesCanopyHeight2018, canopy_imgStack$CaplesCanopyHeight2020)

View(wldf_100)

# plot 450 veg 106

wldf_100$point_d=="450"

wldf_100$point_d[wldf_100$point_d=="754"]

dimnames(wldf_100)[[1]][wldf_100$point_d=="754"]
# [1] "67"

dimnames(wldf_100)[[1]][wldf_100$point_d=="450"]

rowidno <-52

par(mfrow=c(1,1))
par(pty="s")
plot(test_ch18[[rowidno]][,1], test_ch20[[rowidno]][,1], xlab = "Canopy height 2018", ylab = "Canopy height 2020")
abline(0,1)
abline(h = 22, v = 22)
title(main=c("Point 450 canopy height 2018 vs 2020"))
abline(lm(test_ch20[[rowidno]][,1]~test_ch18[[rowidno]][,1]), col="red")

#abline(lm(test_ch20[[rowidno]][,1] ~ 1 + offset(1* test_ch18[[rowidno]][,1])), col="red")

test_ch18[[rowidno]][,1] >= 22 -> ch2018g22

sum(ch2018g22 * 1)

test_ch20[[rowidno]][,1] >= 22 -> ch2022g22

sum(ch2022g22 * 1)

#what about CC
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2018, wldf_100) -> test_cc18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyCover2020, wldf_100) -> test_cc20 # %>% head()

# plot(canopy_imgStack$CaplesCanopyHeight2018, canopy_imgStack$CaplesCanopyHeight2020)

View(wldf_100)

# plot 450 veg 106

wldf_100$point_d=="450"

wldf_100$point_d[wldf_100$point_d=="450"]

dimnames(wldf_100)[[1]][wldf_100$point_d=="450"]
# [1] 52

rowidno <- 52

par(mfrow=c(1,1))
par(mar=c(4,2,3,0))
par(oma=c(5,4,5,0))
par(pty="s")
plot(test_cc18[[rowidno]][,1], test_cc20[[rowidno]][,1], xlab = "Canopy Cover 2018", ylab = "Canopy Cover 2020")
abline(0,1)
# abline(h = 22, v = 22)
title(main=c("Point 450 canopy cover 2018 vs 2020"))
abline(lm(test_cc20[[rowidno]][,1]~test_cc18[[rowidno]][,1]), col="red")

#abline(lm(test_cc20[[rowidno]][,1] ~ 1 + offset(1* test_cc18[[rowidno]][,1])), col="red")

test_cc18[[rowidno]][,1] >= 22 -> ch2018g22

sum(ch2018g22 * 1)

test_cc20[[rowidno]][,1] >= 22 -> ch2022g22

sum(ch2022g22 * 1)



