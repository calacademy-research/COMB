library(MASS)
#y.aru[i,j] ~ dpois(lam*z[i] + ome) 
#site.prob[i] <- lam*z[i]/(lam*z[i]+ome)
jagsResult
# attach(jagsResult$sims.list)
nIter <- length(jagsResult$sims.list$p11)
Tobs <- numeric(nIter)
Tsim <- numeric(nIter)
ySim <- list(nIter)
n <- rowSums(!is.na(data$y.aru[, 1:24])) # n = number of files between 6 and 10 am
count_aru <- rowSums(data$y.aru, na.rm = TRUE)
for (iter in 1:nIter) {
  Tobs[iter] <- sum((sqrt(count_aru) - sqrt(
    (
      jagsResult$sims.list$z[iter,] * jagsResult$sims.list$lam[iter] + jagsResult$sims.list$ome[iter]
    )* n)) ^ 2)
  ySim <-
    rpois(
      84,
      jagsResult$sims.list$z[iter,] * jagsResult$sims.list$lam[iter] + jagsResult$sims.list$ome[iter]
    ) * n 
  Tsim[iter] <-
    sum(
      (sqrt(ySim) - sqrt((
        jagsResult$sims.list$z[iter,] * jagsResult$sims.list$lam[iter] + jagsResult$sims.list$ome[iter]
      ) * n
      )) ^ 2)
}
# detach(jagsResult$sims.list)



MASS::eqscplot(
  Tobs,
  Tsim,
  xlim = range(Tobs, Tsim),
  ylim = range(Tobs, Tsim),
  xlab = "Observed data",
  ylab = "Simulated data"
)
abline(0, 1, lwd = 2, col = "red")
mean(Tsim > Tobs) # the P value

hist(Tsim, breaks = 50)
hist(Tobs, breaks = 50)
