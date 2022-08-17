exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2018, wldf_100) -> test_ch18 # %>% head()
exactextractr::exact_extract(canopy_imgStack$CaplesCanopyHeight2020, wldf_100) -> test_ch20 # %>% head()

plot(canopy_imgStack$CaplesCanopyHeight2018, canopy_imgStack$CaplesCanopyHeight2020)

View(wldf_100)

# plot 450 veg 106

wldf_100$Avian_Poin=="450"

wldf_100$Avian_Poin[wldf_100$Avian_Poin=="754"]

dimnames(wldf_100)[[1]][wldf_100$Avian_Poin=="754"]
[1] "67"

dimnames(wldf_100)[[1]][wldf_100$Avian_Poin=="450"]

rowidno <- 85

par(mfrow=c(1,1))
par(pty="s")
plot(test_ch18[[rowidno]][,1], test_ch20[[rowidno]][,1], xlab = "CH 2018", ylab = "CH 2020")
abline(0,1)
abline(h = 22, v = 22)

abline(lm(test_ch20[[rowidno]][,1]~test_ch18[[rowidno]][,1]), col="red")

#abline(lm(test_ch20[[rowidno]][,1] ~ 1 + offset(1* test_ch18[[rowidno]][,1])), col="red")

test_ch18[[rowidno]][,1] >= 22 -> ch2018g22

sum(ch2018g22 * 1)

test_ch20[[rowidno]][,1] >= 22 -> ch2022g22

sum(ch2022g22 * 1)


