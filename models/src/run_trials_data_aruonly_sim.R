#' run_trials_data: Parallel estimation using multiple hyperparameter settings.
#' This derivative of `run_trials` is meant to run the model on a single speciest with varying
#' the amount of data it is fed (# of PC days/# of ARU recordings)
#'
#' Example usage:
#'
#' ```
#' source("models/src/max_aruonly_arubin_2covars_sim.R")
#' source("models/src/run_trials_data_aruonly_sim.R")
#' startTrials(perTrialHparams("GCKI"), ModelTrial)
#'
#' # Wait for some time
#'
#' getCurrentResults()
#' ```
#'
#' The canonical example is to run JAGS for each of several species in,
#' reducing time to all-species result by a factor of (cores / chains). The
#' code uses a generic list as the hyperparameter data structure.
#'
#' A run of the model under particular hyperparameters should be encapsulated
#' by the caller as a function "trialFn" that takes a single argument, a generic
#' list of parameters.
#'
#' Appendix: Why not purrr?
#'
#' The implementation is parallel and asynchronous, using the promises package.
#' Although this requires polling global state for the results, there is a long
#' tail in the distribution of JAGS time to completion over all possible
#' hyperparameter settings (stragglers). A synchronous API like purrr blocks
#' until all of the trials have completed, which would leave most cores
#' unutilized during the long wait time for the worst straggler to complete.
NULL

library(future)
library(promises)
library(tidyverse)

assumedNumChains <- 2
# Too scary
# plan(multisession, workers = availableCores() / assumedNumChains)

plan(multisession, workers = 48)

#' Generates a hyperparameter sweep over species codes.
#'
#' @return list [list(speciesCode=c, year=2021)] for each species code c.
perTrialHparams <- function(speciesCode) {
  p11 <- 0.5
  p_aru11 <- seq(0, 0.3, by = 0.1)
  p_aru01 <- 0.025 # we should probably not have this parameter...
  psi = seq(0.1, 0.9, by = 0.1)
  
  combs <- expand_grid(0, 1:24, p11, p_aru11, p_aru01, psi) %>% 
    rename(nPC = 1, nARU = 2, p11 = 3, p_aru11 = 4, p_aru01 = 5, psi = 6) %>% 
    as.list()
  
  pmap(combs, function (nARU, nPC, p11, p_aru11, p_aru01, psi) {
    list(speciesCode = speciesCode, year = 2021, nPC = nPC, nARU = nARU, 
         p11 = p11, p_aru11 = p_aru11, p_aru01 = p_aru01, psi = psi)
  }
  )
  
}

#' Global list which will be populated by per-trial calls to a "promise
#' fulfillment" callback.
#'
#' The type is a collection of (params, jagsResult) structures. Formal example:
#'
#' c(list(params={p1}, jagsResult={jr1}), list(params={p2}, jagsResult={jr2}),
#' ...)
trialResults <- list()

#' Starts running a JAGS run function for each setting of hyperparameters.
#'
#' @param hparamsCollection list of lists, where each element specifies a single
#'   trial and is therefore passed as the single argument of runTrial.
#' @param runTrial function(params) -> jagsResult
#'
#' @return NULL
startTrials <- function(hparamsCollection, runTrial) {
  for (hparams in hparamsCollection) {
    promise <- future_promise(
      {
        jagsResult <- runTrial(hparams)
        kv <- list(result = jagsResult)
        # TODO(matt.har.vey): Created a proper (params, jagsResult) data
        # structure and use it to pass param values through to fields of the
        # results frame. The current code is cheating and using a list to pass
        # the species code, but that leaves no room for other params.
        names(kv) <- paste(hparams$nPC, hparams$nARU, hparams$psi, hparams$p11, 
                           hparams$p_aru11, hparams$p_aru01, sep = "_")
        return(kv)
      },
      seed = 123
    )
    then(
      promise,
      onFulfilled = function(kv) {
        trialResults <<- c(trialResults, kv)
      },
      onRejected = function(err) {
        warning(hparams$speciesCode, " failed: ", err)
      }
    )
  }
  NULL
}

#' Generates the "p" value for the pp_checks (posterior predictive checks)
#' Given the name of the simulated values and the observed values
#' 
#' @param jagsResults as returned from a call to jagsUI::jags()
#' @param obs_sim is a named vector: c("Dobs" = "Dsim", "Tobs" = "Tsim")
#' 
#' @return list(Dobs_Dsim={}, ...)
pp_check <- function(jagsResult, obs_sim){
  simple_check <- function(observed, simulated){
    obs <- c(mcmc_to_mat(jagsResult$samples, observed))
    sim <- c(mcmc_to_mat(jagsResult$samples, simulated))
    bpval <- mean(sim > obs)
  }
  pp_vals <- map2(
    names(obs_sim), obs_sim,
    simple_check
  )
  names(pp_vals) <- paste(names(obs_sim), obs_sim, sep = "_")
  unlist(pp_vals)
}

#' Finds some sort of statistic (eg. variance) for the posterior samples
#' for a parameter (or multiple parameters)
#' 
#' @param jagsResult as returned from a call to jagsUI::jags()
#' @param FUN Functions to be used on the posterior samples
#' @param parameter Vector of parameters to have functions preformed on
#' 
#' @return list()
post_stats <- function(jagsResult, FUN, parameter){
  x <- map(parameter, 
           function(.){jagsResult$sims.list[[.]] %>% FUN}
  ) %>% unlist()
  
  names(x) <- paste0(parameter, "_var") # change 'var' to other fun later
  x
}

#' Extracts a list of posterior parameter means and Rhats.
#'
#' @param jagsResult as returned from a call to jagsUI::jags()
#'
#' @return list(psi={}, ..., rhat_psi={}, ...)
collectEstimates <- function(jagsResult) {
  means <- unlist(jagsResult$mean)
  rhats <- unlist(jagsResult$Rhat)
  names(rhats) <- as.character(
    map(names(rhats), function(n) {
      paste("rhat", n, sep = "_")
    })
  )
  pp <- pp_check(jagsResult, c("Dobs" = "Dsim", 
                               "TvegObs" = "TvegSim", 
                               "T_pc_obs" = "T_pc_sim",
                               "T_aru_obs" = "T_aru_sim"))
  stats <- post_stats(jagsResult, FUN = function(.){1/var(.)}, "mean_psi")
  append(stats, means) %>% append(pp) %>% append(rhats)
}



#' Returns a tibble of (params, estimates) for trials that have finished.
#'
#' The design is that this function be called manually and periodically to
#' inspect results of the trials that have completed so far.
#'
#' @return tibble where each row collects hyperparameter values and posterior
#'   means of the model parameters.
getCurrentResults <- function() {
  table <- t(sapply(trialResults, collectEstimates))
  data.frame(table) %>% rownames_to_column("params") %>% 
    separate(params, c("nPC", "nARU", "psi", "p11", "p_aru11", "p_aru01")) %>% 
    mutate(nPC = as.numeric(nPC), nARU = as.numeric(nARU), 
           psi = as.numeric(psi), p11 = as.numeric(p11), 
           p_aru11 = as.numeric(p_aru11), p_aru01 = as.numeric(p_aru01)
    )
}

# x <- getCurrentResults()

## Plotting the two Gaussian Distributions for the Scores
# ggplot(x, aes(x = seq(-4, 5, length.out = nrow(x)))) + 
#   geom_area(stat = "function", fun = dnorm, args = c(mean = x$mu1[1], sd = (1/sqrt(x$sigma[1]))), fill = "blue", alpha = 0.5) + 
#   geom_area(stat = "function", fun = dnorm, args = c(mean = x$mu2[1], sd = (1/sqrt(x$sigma[1]))), fill = "red", alpha = 0.5) + 
#   labs(title = "The Two Gaussian Distributions for the Gaussian Mixture Model") + 
#   ylab("density") +xlab("logit")
