
# attach(jagsResult_WEWP$sims.list)
nIter <- length(jagsResult_WEWP$sims.list$p11)
Tobs <- numeric(nIter)
Tsim <- numeric(nIter)
n <- rowSums(!is.na(data$y.ind[, 1:4])) # n = number of visits
y <- rowSums(data$y.ind[, 1:4], na.rm = TRUE)
for (iter in 1:nIter) {
  Tobs[iter] <- sum((sqrt(y) - sqrt((jagsResult_WEWP$sims.list$z[iter, ] * jagsResult_WEWP$sims.list$p11[iter] + (1 - jagsResult_WEWP$sims.list$z[iter, ]) * jagsResult_WEWP$sims.list$p10[iter]) * n))^2)
  ySim <- rbinom(84, n, jagsResult_WEWP$sims.list$z[iter, ] * jagsResult_WEWP$sims.list$p11[iter] + (1 - jagsResult_WEWP$sims.list$z[iter, ]) * jagsResult_WEWP$sims.list$p10[iter])
  Tsim[iter] <- sum((sqrt(ySim) - sqrt((jagsResult_WEWP$sims.list$z[iter, ] * jagsResult_WEWP$sims.list$p11[iter] + (1 - jagsResult_WEWP$sims.list$z[iter, ]) * jagsResult_WEWP$sims.list$p10[iter]) * n))^2)
}
# detach(jagsResult_WEWP$sims.list)

MASS::eqscplot(Tobs, Tsim,
               xlim = range(Tobs, Tsim), ylim = range(Tobs, Tsim),
               xlab = "Observed data", ylab = "Simulated data"
)
abline(0, 1, lwd = 2, col = "red")
mean(Tsim > Tobs) # the P value

hist(Tsim, breaks = 50)
hist(Tobs, breaks = 50)
