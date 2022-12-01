#RAVG_cbi4_vs_cbi.R
#ravg comparisons cbi4 vs cbi (old vs new)
par(mfrow=c(1,2))

cbi4 <- RAVG_Caples_imgStack$ca3872412014620191010_20181118_20191118_rdnbr_cbi4
cbi <- RAVG_Caples_imgStack$ca3872412014620191010_20181118_20191118_rdnbr_cbi

cbi4 %>% 
  cut(.,c(-1,0,1,3,4)) -> cbi4_cuts

cbi %>% 
  cut(.,c(0,0.1,1.25,2.25,3.0)) -> cbi_cuts

hist(cbi4_cuts, main="CBI4")
hist(cbi_cuts, main="CBI")

cbis <- raster::stack(cbi, cbi4)

plot(cbis, maxpixels=1e5)



