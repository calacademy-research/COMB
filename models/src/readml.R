#' readml library
#'
#' This file defines functions for reading machine learning model outputs in
#' the format used by COMB and for aggregating and structuring them into the
#' type of data lists used by `AHMbook`.
NULL

library(lubridate)
library(tibble)

#' Reads and structures machine learning model outputs
#'
#' Provides a one-line way to get a data list acceptable to `jags`.
#'
#' The main output is a multidimensional array of counts (y.aru) with dimensions
#' being one of:
#'
#'   - [species, years, sites, visits]
#'   - [species, sites, visits]
#'   - [years, sites, visits]
#'   - [sites, visits]
#'
#' according to whether the squeeze param is true and to the lengths of the
#' species and years params (Note the convention that kept axes are in the same
#' order as if all axes had been present.)
#'
#' @param species Vector of species codes in the same format as the input CSV.
#'   Order determines the indexing of the species axis.
#' @param years Vector of integer years. Order determines the indexing of the
#'   year axis.
#' @param beginTime Optional `lubridate::duration` offset since midnight for
#'   which all earlier-in-day scores should be dropped.
#' @param endTime Optional `lubridate::duration` offset since midnight for which
#'   all later-in-day scores should be dropped.
#' @param visitAggregation 'file' or 'day', indicating which scope to group at
#'   before sorting the groups by time onto the visits axis.
#' @param visitLimit Only this many visits will be kept, per the definition and
#'   ordering determined by visitAggregation.
#' @param thresholdOptions list with two names
#'     - value: numeric scalar to use the same threshold for all species or
#'       vector of per-species thresholds, indexed the same as the species
#'       param.
#'     - is.quantile: Whether the threshold values should be interpreted as
#'       quantiles of all scores.
#' @param squeeze Whether to drop the year or species axes when their size is 1.
#'   This lets this same function work for any rank of counts matrix assumed by
#'   the downstream JAGS model.
#'
#' @return list of data, including all variable names used in JAGS models
#'   - Dimensions
#'     - nsites: The number of sites. (Length of the site axis in y.aru.)
#'     - nsurveys.aru: The number of ARU visits.
#'   - Counts
#'     - y.aru: Multidimensional array of above-threshold scores.
#'   - Scores
#'     - score: Vector of scores, the number of which tabulates to the counts in
#'       y.aru. (Can be viewed as the values list of a sparse matrix.)
#'     - siteid: Site indices corresponding to score. (Can be viewed as row
#'       indices into a sparse matrix.)
#'     - occid: ARU "visit" index corresponding to score. (Can be viewed as
#'       column indices into a sparse matrix.)
#'     - nsamples: length(scores), convenience alias for JAGS models.
#'   - Indexing (for reference)
#'     - point_numbers: Vector of point numbers, corresponding to the point axis
#'       of y.aru.
#'     - visit_times: Vector of POSIXct instants, corresponding to the visit
#'       axis of y.aru. (Intended to aid interpretation of visit indices.
#'       The values are nominal, not precise, for example start-of-file or
#'       midnight.)
#'
#' @export
readMLOutputs <- function(species, years, beginTime, endTime,
                          visitAggregation = "file", visitLimit = NA,
                          thresholdOptions, squeeze = T) {
  scores <- readMLTibble(species, years, beginTime, endTime)
  structureMLData(scores, thresholdOptions)
  stop("not implemented")
}

#' Reads long (aka tall) machine learning model outputs
#'
#' Reads from well known "tall" CSV file(s), filters to a specified time-of-day
#' interval, and concatenates to a single tibble.
#'
#' Restricting the read to subsets of species and years given as params enables
#' targeted, and therefore quicker, reads.
#'
#' @param species Vector of species codes to include in the result.
#' @param years Vector of numeric years to include in the result.
#' @param beginTime Optional `lubridate::duration` offset since midnight for
#'   which all earlier-in-day scores should be dropped.
#' @param endTime Optional `lubridate::duration` offset since midnight for which
#'   all later-in-day scores should be dropped.
#'
#' @return tibble with columns [Date_Time, Point, Start_Time, Species, Logit]
#'
#' @export
readMLTibble <- function(species, years, beginTime, endTime) {
  stop("not implemented")
}

#' Structures machine learning model outputs as JAGS data
#'
#' Adapts the output of readMLTibble to the `data` for `jags`.
#'
#' @param scores tibble with columns [Date_Time, Point, Start_Time, Species,
#'   Logit] (as returned by readMLTibble).
#' @param visitAggregation See corresponding parameter of readMLOutputs.
#' @param visitLimit See corresponding parameter of readMLOutputs.
#' @param thresholdOptions See corresponding parameter of readMLOutputs.
#' @param squeeze See corresponding parameter of readMLOutputs.
#'
#' @return list of `data` for `jags`. See the return of readMLOutputs.
#'
#' @export
structureMLData <- function(scores, visitAggregation = "file", visitLimit = NA,
                            thresholdOptions, squeeze = T) {
  stop("not implemented")
}
