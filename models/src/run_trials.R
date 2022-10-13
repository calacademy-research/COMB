#' run_trials: Parallel estimation using multiple hyperparameter settings.
#'
#' Example usage:
#'
#' ```
#' source('models/src/combined.R')
#' source('models/src/run_trials.R')
#' startTrials(perSpeciesHparams(), singleSpeciesCombined)
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
#' by the caller as a funtion "trialFn" that takes a single argument, a generic
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

plan(multisession, workers = 32)

#' Generates a hyperparameter sweep over species codes.
#'
#' @return list [list(speciesCode=c, year=2021)] for each species code c.
perSpeciesHparams <- function() {
  guilds <- read_csv(
    "models/input/bird_guilds.csv",
    col_names = c("name", "code6", "code4", "guild"),
    col_types = cols(
      name = col_character(),
      code6 = col_character(),
      code4 = col_character(),
      guild = col_character()
    )
  )
  map(guilds$code4, function(c) {
    list(speciesCode = c, year = 2021, nARU = 24, nPC = 3)
  })
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
        names(kv) <- hparams$speciesCode
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
  append(means, rhats)
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
  data.frame(table) %>% rownames_to_column("species")
}
