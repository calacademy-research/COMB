## run_simulation.R
## Script to run simulation and aggregate the outputs

# Import script that contains the function to run the model
source("models/src/max_pc_scores_arubin_2covars_sim.R")
library(future)
library(promises)

plan(multisession, workers = 2)

simResults <- list()

params <- list(
  nARU = 24,
  nPC = 3, 
  p11 = 0.5, 
  p_aru11 = 0.4, 
  p_aru01 = 0.025, 
  beta0 = 0, 
  beta1 = 0
)

#' run_sim
#' @param params list of params to be fed into the simulated data set
#' @params n_iter integer number of iterations the model should be run
run_sim <- function(params, n_iter) {
  for (i in 1:n_iter) {
    promise <- future_promise({
      result <- ModelTrial(params)
      kv <- list(result)
      return(kv)
    }, seed = T)
    then(
      promise,
      onFulfilled = function(kv) {
        simResults <<- c(simResults, kv)
      },
      onRejected = function(err) {
        warning("trial failed", err)
      }
    )
  }
}

# means <- map(1:length(simResults), 
             ~ simResults[[.x]][["mean"]][["beta0"]]
             ) %>% unlist
