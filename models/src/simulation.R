## simulation.R
## Script to simulate data to 'test' the model

library(tidyverse)
library(truncnorm)
source("comb_functions.R")

#' simulation
#' @param species Species Code that we will simulate data for
#' @param psi Theoretical probability of occupancy for the given species
#' @param n_visits Number of visits by human point counters
#' @param n_points Number of points
#' @param n_recordings Number of recordings at each point by an ARU
#' @param mu Vector containing the means for the two Gaussian distributions of
#' the scores. The first is the negative mu and the second is the positive.
#' @param sigma Vector containing the standard deviations of the two Gaussian
#' distributions of the scores. The first is the negative mu and the second is
#' the positive.
#' @return List similar to that of `model_read_lib` that JAGS can read

simulation <- function(species, psi, p11, p_aru11, p_aru01, n_visits = 3,
                       n_points = 82, n_recordings = 24, mu, sigma) {
  ## Point Count Simulation ----------------------------------------------------
  # Generate some Bernoulli trials
  set.seed(123)

  z <- rbinom(n_points, 1, psi) # occupancy states

  # p <- psi * p11 # detection probability

  y.ind <- matrix(NA, nrow = n_points, ncol = n_visits)

  for (i in 1:n_points) {
    pr_yequals1 <- p11 * z[i]
    y.ind[i, ] <- rbinom(n_visits, 1, pr_yequals1)
  }

  # rbinom(n_visits * n_points, size = 1, prob = p) %>%
  # matrix(nrow = n_points, ncol = n_visits)

  ## ARU Simulation ------------------------------------------------------------

  # We will first define the `y.aru` which is a matrix similar to y.ind
  # and has binary data on if the ARU has a detection of the species above an
  # arbitrary logit
  nsamples <- n_recordings * n_points

  y.aru <- matrix(NA, nrow = n_points, ncol = n_recordings)

  for (i in 1:n_points) {
    p_aru <- z[i] * p_aru11 + p_aru01
    y.aru[i, ] <- rbinom(n_recordings, 1, p_aru)
  }


  # Next we will simulate the scores

  # The scores are distributed as two truncated Gaussian distributions
  # with parameters mu and sigma.

  scores <- c()

  # for each of the positive pts, simulate a vector of scores with len=n_recordings
  for (i in 1:n_points) {
    if (n_recordings != 0) {
      for (j in 1:n_recordings) {
        if (z[i] == 1) {
          scores[(i - 1) * n_recordings + j] <- rnorm(1, mean = mu[2], sd = sigma[2])
        }
        if (z[i] == 0) {
          scores[(i - 1) * n_recordings + j] <- rnorm(1, mean = mu[1], sd = sigma[1])
        }
      }
    }
  }

  sites <- rep(1:n_points, each = n_recordings)


  # Combining for JAGS ---------------------------------------------------------

  data <- list(
    "nsites" = n_points,
    "nsurveys.pc" = n_visits,
    "nsurveys.aru" = n_recordings,
    "y.ind" = y.ind,
    "y.aru" = y.aru,
    "nsamples" = length(scores),
    "siteid" = sites,
    "score" = scores
  )
  data
}

# data <- simulation("GCKI", n_points = 80, psi = 0.62, p11 = 0.6, p_aru11 = 0.7,
#                    p_aru01 = 0.05, mu = c(-2, 1.5), sigma = c(0.8, 2.3),
#                    n_visits = 3, n_recordings = 24)
