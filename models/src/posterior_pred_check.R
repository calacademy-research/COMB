library(MASS)

jagsResult
# attach(jagsResult$sims.list)
nIter <- length(jagsResult$sims.list$p11)
Tobs <- numeric(nIter)
Tsim <- numeric(nIter)
n <- rowSums(!is.na(data$y.ind[, 1:3])) # n = number of visits
y <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
for (iter in 1:nIter) {
  Tobs[iter] <- sum((sqrt(y) - sqrt((jagsResult$sims.list$z[iter, ] * jagsResult$sims.list$p11[iter] + (1 - jagsResult$sims.list$z[iter, ]) * jagsResult$sims.list$p10[iter]) * n))^2)
  ySim <- rbinom(84, n, jagsResult$sims.list$z[iter, ] * jagsResult$sims.list$p11[iter] + (1 - jagsResult$sims.list$z[iter, ]) * jagsResult$sims.list$p10[iter])
  Tsim[iter] <- sum((sqrt(ySim) - sqrt((jagsResult$sims.list$z[iter, ] * jagsResult$sims.list$p11[iter] + (1 - jagsResult$sims.list$z[iter, ]) * jagsResult$sims.list$p10[iter]) * n))^2)
}
# detach(jagsResult$sims.list)

MASS::eqscplot(Tobs, Tsim,
               xlim = range(Tobs, Tsim), ylim = range(Tobs, Tsim),
               xlab = "Observed data", ylab = "Simulated data"
)
abline(0, 1, lwd = 2, col = "red")
mean(Tsim > Tobs) # the P value

hist(Tsim, breaks = 50)
hist(Tobs, breaks = 50)
