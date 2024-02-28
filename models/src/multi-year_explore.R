# multi-year explore

# BBWO 

bsum <-data.frame(jagsResult_BBWO_multiyear$summary)
bsum$measure <- rownames(bsum)


beta0 <- bsum[2:5,]
beta0$year <- c(2018,2019, 2020, 2021)
beta0$prepost <- c("pre", "pre", "post", "post")

ggplot(beta0) + 
  geom_point(aes(x=year, y=plogis(mean), color=prepost)) +
  geom_errorbar(aes(x=year, ymin=plogis(X2.5.+), ymax=plogis(X97.5.), color=prepost))












# scores [under construction] ---------------------------------------------


scorepre <- jagsData$score[jagsData$yearid==1|jagsData$yearid==2]
score19 <- jagsData$score[jagsData$yearid==2]
scorepost <- jagsData$score[jagsData$yearid==3|jagsData$yearid==4]
score21 <- jagsData$score[jagsData$yearid==4]
# score distributions by year
# scoreyr <- 
# for (i in 1:jagsData$nyears) {
#   scoreyr$i <- jagsData$score[jagsData$yearid==i]
# }

pre <- hist(scorepre[scorepre>0], breaks=100, xlab = "ML score distribution", xlim = c(-3,15))
post <- hist(scorepost[scorepost>0], breaks=100, xlab = "ML score distribution", xlim = c(-3,15))
plot(pre, col=alpha("pink", 0.5), add=T)
plot(post, col=alpha("forestgreen",  0.5))
