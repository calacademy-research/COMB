#' model_read_lib: Functions to read observations used in JAGS models
#'
#' readCombined is the "do everything" function.
#'
#' This file defines functions for reading scores indexed by some subset of
#' (Species, Year, Point, Visit). Some functions deal with "scores" generically.
#' Higher-level functions for reading point count data, ML outputs, and the
#' combination of the two, call into the generic functions and collect the
#' output under qualified names.
#'
#' The high-level functions aim to make their outputs data lists, where the
#' names are a superset of those used in JAGS models from `AHMbook`.
NULL

library(fst)
library(dplyr)
library(here)
library(lubridate)
library(readr)
library(stringr)
library(tibble)
library(data.table)

# All of these are symlinks, which adds a layer of indirection, so we don't do
# any drive_sync here.
latlongPath <- here("models/input/latlong.csv")
aru2pointPath <- here("models/input/aru2point.csv")
dataMlPath <- here("models/input/dataML_tall.csv")
pointCountsPath <- here("models/input/PC_delinted.csv")

# The functions in this library are ordered approximately top-down, from those
# most likely to be called by a dependent script to those most purely internal /
# supportive of other functions in this library.

#' Reads point counts and ML model outputs and structures for JAGS
#'
#' The output is a data list, the most notable names of which are:
#'
#'   - `y.ind`: binary point counts
#'   - `y.pc`: raw point counts
#'   - `y.aru`: ML detection counts
#'
#' The values of there are all multidimensional arrays of the same shape, one
#' of:
#'
#'   - [species, years, sites, visits]
#'   - [species, sites, visits]
#'   - [years, sites, visits]
#'   - [sites, visits]
#'
#' according to whether the squeeze param is TRUE and to the lengths of the
#' species and years params (Note that the squeezed versions keep the same
#' ordering of dimensions as the full version.)
#'
#' The returned data list also includes many other names conventional in the
#' JAGS models in `AHMbook`.
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
#'   before sorting the groups by time onto the visits axis. (Only applies to
#'   ARU visits.)
#' @param visitLimit Only this many visits will be kept, per the definition and
#'   ordering determined by visitAggregation. (Only applies to ARU visits.)
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
#' @return list of data, including all variable names used in JAGS models as
#'   well as additional values for reference.
#'   - Indices
#'     - indices: list with names (species, year, point, visit.pc, visit.aru),
#'       where the values are data.tables that give the 1:1 mapping from indices
#'       to the corresponding values ((Species, Species_Index), etc.)
#'   - Dimensions
#'     - nspecies: Maximum index on the species axis.
#'     - nyears: Maximum index on the year axis.
#'     - nsites: Maximum index on the point axis.
#'     - nsurveys.pc: The number of point counts. (Inner dimension of y.pc and
#'       y.ind.)
#'     - nsurveys.aru: The number ARU "visits" count. (Inner dimension of y.aru.
#'       See also visitAggregation.)
#'   - Counts
#'     - y.ind: Binarized point counts (occupancy per observers).
#'     - y.pc: Raw point counts.
#'     - y.aru: Count of above-threshold ML scores.
#'   - Scores
#'     - nsamples: The number of above-threshold ML scores.
#'     - speciesid: Species indices corresponding to score.
#'     - yearid: Year indices corresponding to score.
#'     - siteid: Site indices corresponding to score.
#'     - occid: ARU "visit" index corresponding to score.
#'     - score: Vector of scores. (As a check, (len(score) == sum(y.aru)).)
#'
#' @export
readCombined <- function(species, years, beginTime = NA, endTime = dhours(10),
                         visitAggregation = "file", visitLimit = NA, PCvisitlimit = NA, 
                         thresholdOptions = list(
                           value = -2.0,
                           is.quantile = F
                         ),
                         squeeze = T) {
  outerIndices <- buildOuterIndices(species, years)
  pointCountData <- readPointCounts(outerIndices, PCvisitlimit = PCvisitlimit, squeeze = squeeze)
  aruData <- readML(
    outerIndices,
    beginTime = beginTime, endTime = endTime,
    visitAggregation = visitAggregation, visitLimit = visitLimit,
    thresholdOptions = thresholdOptions, squeeze = squeeze
  )
  combineJagsData(pointCountData, aruData)
}

#' Combines point counts with ML outputs based on ARU data
#'
#' This is a restructuring method that combines the data lists returned by
#' point-count-specific and ML-specific functions. Mainly, this means checking
#' the consistency of the outer indices and qualifying names.
#'
#' @param pointCountData JAGS data list based on point counts. Usually the
#'   output of readPointCounts.
#' @param aruData JAGS data list based on ML detections from ARUs. Usually the
#'   output of readML.
#'
#' @return list of data for JAGS. It is roughly a union of pointCountData and
#'   aruData, with common names keps only once and with qualifying suffixes
#'   appended to names that collide (y, nvisits). For full detail, see the
#'   documentation of readComined.
#'
#' @export
combineJagsData <- function(pointCountData, aruData) {
  # Make sure that point count and ML indices correspond 1:1, in order.
  outerIndexNames <- c("species", "year", "point")
  getOuterIndices <- function(data) {
    data$indices[outerIndexNames]
  }
  outerIndices <- getOuterIndices(pointCountData)
  stopifnot(identical(outerIndices, getOuterIndices(aruData)))
  
  s <- aruData$sparseScore # alias for readability
  
  list(
    indices = c(
      outerIndices,
      list(
        visit.pc = pointCountData$indices$visit,
        visit.aru = aruData$indices$visit
      )
    ),
    nspecies = max(outerIndices$species$Species_Index),
    nyears = max(outerIndices$year$Year_Index),
    nsites = max(outerIndices$point$Point_Index),
    nsurveys.pc = max(pointCountData$indices$visit$Visit_Index),
    nsurveys.aru = max(aruData$indices$visit$Visit_Index),
    y.ind = as.matrix(pointCountData$y),
    y.pc = pointCountData$y.raw,
    y.aru = as.matrix(aruData$y),
    # ARU scores (sparse)
    nsamples = nrow(s),
    speciesid = s$Species_Index,
    yearid = s$Year_Index,
    siteid = s$Point_Index,
    occid = s$Visit_Index,
    score = s$Score
  )
}

#' Reads point counts and builds a JAGS data list
#'
#' This is the point-count-specific branch of readCombined. See the
#' documentation there for full details.
#'
#' @param outerIndices List of the form returned by buildOuterIndices. Passing
#'   the same to readML ensures consistent intrepretation of [species, year,
#'   point] indices.
#'
#' @return return JAGS data list as returned by readCombined except that y.ind
#'   is named just y, nsurveys.pc is named just nsurveys, and ARU-specific names
#'   are not present.
#'
#' @export
readPointCounts <- function(outerIndices, PCvisitlimit, squeeze = T) {
  counts <- read_csv(
    pointCountsPath,
    col_types = cols(
      observer_fk = col_character(),
      birdCode_fk = col_character(),
      abun = col_integer(),
      point_ID_fk = col_integer(),
      year = col_integer(),
      visit = col_integer(),
      DateTime = col_datetime()
    )
  ) %>% mutate(
    Species = birdCode_fk, Year = year(DateTime),
    Point = point_ID_fk, Visit = visit, Score = abun
  ) %>% filter(Visit <= PCvisitlimit)
  # arrange(Point, Year, Visit) %>% 
  # group_by(Point, Year) %>% 
  # mutate(
  #  Visit = seq_along(Visit)
  # ) #%>% 
  # filter(Visit <= PCvisitlimit) %>% ungroup()
  
  # (y == 1 if occupied) can be had by treating counts as "scores" and setting
  # the threshold to 0.
  visits <- counts %>%
    select(Year, Point, Visit) %>%
    distinct()
  scores <- counts %>%
    filter(Score > 0) %>%
    select(Species, Year, Point, Visit, Score)
  countsData <- structureForJags(outerIndices, visits, scores,
                                 visitLimit = NA,
                                 squeeze = squeeze
  )
  
  # Add raw counts to the data list, to give the option of modeling a rate of
  # observer counts.
  indices <- buildFullIndices(outerIndices, visits, visitLimit = NA)
  sparseRawCounts <- counts %>%
    inner_join(indices$full, by = c("Species", "Year", "Point", "Visit")) %>%
    select(Species_Index, Year_Index, Point_Index, Score)
  y.raw <- sparseToDense(sparseRawCounts, indices$full)
  dimnames(y.raw)[[1]] <- indices$species$Species
  dimnames(y.raw)[[2]] <- indices$year$Year
  dimnames(y.raw)[[3]] <- indices$point$Point
  
  c(countsData, list(y.raw = y.raw))
}

#' Reads and structures machine learning model outputs
#'
#' This is the ARU-specific branch of readCombined. See the documentation there
#' for full details.
#'
#' @param outerIndices List of the form returned by buildOuterIndices. Passing
#'   the same to readPointCounts ensures consistent intrepretation of [species,
#'   year, point] indices.
#'
#' @return return JAGS data list as returned by readCombined except that y.aru
#'   is named just y, nsurveys.aru is named just nsurveys, and point
#'   count-specific names are not present.
#'
#' @export
readML <- function(outerIndices, beginTime = NA, endTime = dhours(10),
                   visitAggregation = "file", visitLimit = NA,
                   thresholdOptions = list(value = -2.0, is.quantile = F),
                   squeeze = T) {
  aru2point <- readAru2point()
  mlTibble <- readDataMl(
    outerIndices$species$Species,
    outerIndices$year$Year,
    beginTime,
    endTime
  )
  
  # Visits
  if (visitAggregation == "file") {
    addVisitKeys <- function(.) {
      mutate(., Year = year(Date_Time), Visit = Date_Time)
    }
  } else if (visitAggregation == "day") {
    addVisitKeys <- function(.) {
      mutate(., Year = year(Date_Time), Visit = yday(Date_Time))
    }
  }
  visits <- aru2point %>%
    mutate(Date_Time = parse_date_time(
      str_extract(filename, "\\d{8}_\\d{6}"), "%Y%m%d_%H%M%S"
    )) %>%
    filterTimeOfDay(beginTime = beginTime, endTime = endTime) %>%
    select(Point = point, Date_Time) %>%
    addVisitKeys() %>%
    select(Year, Point, Visit, Date_Time)
  
  # Scores
  if (length(thresholdOptions$value) != 1) {
    # TODO: Handle per-species array of thresholds.
    stop("per-species thresholds are not implemented")
  }
  if (thresholdOptions$is.quantile) {
    stop("quantile thresholds are not implemented")
  }
  threshold <- thresholdOptions$value
  species <- outerIndices$species$Species
  years <- outerIndices$year$Year
  scores <- mlTibble %>%
    #filter(Score > threshold) %>% #`commented out to keep all scores`
    addVisitKeys() %>%
    select(Species, Year, Point, Visit, Score)
  
  structureForJags(outerIndices, visits, scores,
                   visitLimit = visitLimit,
                   squeeze = squeeze, 
                   threshold = threshold
  )
}

#' Builds a JAGS data list
#'
#' This adds a visits dimension to the indices, build a "y matrix" that
#' tabulates score counts along the index dimensions and gives names that match
#' `AHMbook` examples to the columns of the scores table.
#'
#' @param outerIndices List of value-to-index tables as returned by
#'   buildOuterIndices. Used to keep index meanings consistent.
#' @param visits Table with (Year, Point, Visit) columns that will be sorted and
#'   to assign a Visit_Index scoped to each (Year, Point).
#' @param visitLimit Optional maximum Visit_Index to consider. Visits assigned
#'   higher indices will be silently ignored.
#' @param scores Table of above-threshold scores with columns (Species, Year,
#'   Point, Visit, Score). Counts of these become the "y matrix."
#' @param squeeze Whether to drop the year or species axes when their size is 1.
#'   This lets this same function work for any rank of counts matrix assumed by
#'   the downstream JAGS model.
#'
#' @return JAGS data list as returned by readCombined except that the names are
#'   unqualified.
#'
#' @export
structureForJags <- function(outerIndices, visits, visitLimit = NA, scores,
                             squeeze = T, threshold = 0.5) {
  indices <- buildFullIndices(outerIndices, visits, visitLimit = visitLimit)
  
  # Scores
  sparseScores <- scores %>%
    select(Species, Year, Point, Visit, Score) %>%
    inner_join(indices$full, by = c("Species", "Year", "Point", "Visit")) %>%
    select(Species_Index, Year_Index, Point_Index, Visit_Index, Score)
  
  # Counts
  groupByIndices <- function(.) {
    group_by(., Species_Index, Year_Index, Point_Index, Visit_Index)
  }
  initialCounts <- indices$full %>% # visited points start from 0
    groupByIndices() %>%
    summarise(Count = 0, .groups = "drop")
  scoreCounts <- sparseScores %>%
    groupByIndices() %>%
    filter(Score > threshold) %>% # TODO make this `0.5` relate to the threshold
    summarise(Count = n(), .groups = "drop")
  sparseCounts <- rbind(initialCounts, scoreCounts) %>%
    groupByIndices() %>% 
    summarise(Count = sum(Count), .groups = "drop")
  
  y.full <- sparseToDense(sparseCounts, indices$full)
  dimnames(y.full)[[1]] <- indices$species$Species
  dimnames(y.full)[[2]] <- indices$year$Year
  dimnames(y.full)[[3]] <- indices$point$Point
  
  if (squeeze) {
    y <- drop(y.full)
  } else {
    y <- y.full
  }
  
  list(
    # for reference
    indices = indices,
    sparseScores = sparseScores,
    sparseCounts = sparseCounts,
    
    # for JAGS
    nspecies = dim(y.full)[1],
    nyears = dim(y.full)[2],
    nsites = dim(y.full)[3],
    nsurveys = dim(y.full)[4],
    y = y,
    nsamples = nrow(sparseScores),
    speciesid = sparseScores$Species_Index,
    yearid = sparseScores$Year_Index,
    siteid = sparseScores$Point_Index,
    occid = sparseScores$Visit_Index,
    score = sparseScores$Score
  )
}

#' Creates value-to-index tables for the (species, year, point) axes
#'
#' Unlike the visits axis, these axes are independent of survey method, i.e.
#' whether or how many times a point was surveyed. Therefore the indexing into
#' these axes (ordering of dimnames) can match for ARU and point count
#' matrices. This function builds and returns a data structure specifying that
#' ordering that is passed to both ARU and point count read methods and can also
#' be used for reference.
#'
#' @param species vector of species codes whose order determines species axis
#'   indices.
#' @param years vector of years whose order determines year axis indices.
#'
#' @return list of tables (species = (Species, Species_Index), etc.)
#'
#' @export
buildOuterIndices <- function(species, years) {
  latlong <- readLatlong()
  
  species <- tibble(Species = species)
  speciesIndices <- species %>% mutate(Species_Index = seq_along(Species))
  
  years <- tibble(Year = years)
  yearIndices <- years %>% mutate(Year_Index = seq_along(Year))
  
  points <- latlong %>%
    select(Point) %>%
    distinct() %>%
    arrange(Point)
  pointIndices <- points %>% mutate(Point_Index = seq_along(Point))
  
  list(
    species = speciesIndices,
    year = yearIndices,
    point = pointIndices
  )
}

#' Adds $visit and $full to a list of index tables
#'
#' The length and of the visits axis and the interpretation of its indices can
#' vary with survey method. This funtion takes a table of visits, assigns them
#' visit indices, and returns an extended version of the indices data structure
#' given as a param.
#'
#' @param outerIndices list of tables (species = (Species, Species_Index),
#'   year = (Year, Year_Index), point = (Point, Point_Index)) specifying the
#'   value-to-index assignments for the axes that are independent of survey
#'   method.
#' @param visits table with columns (Year, Point, Visit) for a particular survey
#'   method.
#' @param visitLimit Optional. When specified, visits assigned greater indices
#'   will be omitted.
#'
#' @return list like the outerIndices param but extended by (visit = (Visit,
#'   Visit_Index), full = join of all other indices). Any table with all of
#'   (Species, Year, Point, Visit) can be joined on those values to indices$full
#'   in order to get the corresponding (Species_Index, Year_Index, Point_Index,
#'   Visit_Index).
#'
#' @export
buildFullIndices <- function(outerIndices, visits, visitLimit = NA) {
  visits <- visits %>%
    select(Year, Point, Visit) %>%
    distinct() %>%
    inner_join(outerIndices$year, by = "Year")
  visitIndices <- visits %>%
    group_by(Year, Point) %>%
    arrange(Year, Point, Visit) %>% # order important for seq_along
    mutate(Visit_Index = seq_along(Visit))
  if (!is.na(visitLimit)) {
    visitIndices <- visitIndices %>% filter(Visit_Index <= visitLimit)
  }
  
  fullIndices <- visitIndices %>%
    inner_join(outerIndices$point, by = "Point") %>%
    full_join(outerIndices$species, by = character()) %>%
    select(
      Species_Index, Year_Index, Point_Index, Visit_Index,
      Species, Year, Point, Visit
    ) %>%
    arrange(Species_Index, Year_Index, Point_Index, Visit_Index)
  
  c(
    outerIndices,
    list(
      visit = visitIndices,
      full = fullIndices
    )
  )
}

#' Reads machine learning model outputs
#'
#' Reads from "tall" score CSV file(s), filters to a specified time-of-day
#' interval, and concatenates if necessary.
#'
#' @param species Vector of species codes to include in the result.
#' @param years Vector of numeric years to include in the result.
#' @param beginTime Optional `lubridate::duration` offset since midnight for
#'   which all earlier-in-day scores should be dropped.
#' @param endTime Optional `lubridate::duration` offset since midnight for which
#'   all later-in-day scores should be dropped.
#'
#' @return a single tibble with columns [Species, Point, Date_Time, Score]
#'
#' @export
readDataMl <- function(species, years, beginTime = NA, endTime = dhours(10)) {
  fread(
    dataMlPath
  ) %>%
    select(Species = species, Point = point, Date_Time, Score = logit) %>%
    filter(Species %in% species) %>%
    filter(year(Date_Time) %in% years) %>%
    filterTimeOfDay(beginTime = beginTime, endTime = endTime)
}

#' Reads latlong.csv
#'
#' This is the authoritative manifest of points. It provides a consistent
#' assignment of Point_Index for all arrays indexed by point.
#'
#' @return return tibble of (point, latitude, longitude), in the original order
#'   of latlong.csv.
#'
#' @export
readLatlong <- function() {
  read_csv(
    latlongPath,
    col_names = c("Point", "Latitude", "Longitude"),
    col_types = list(
      Point = col_integer(),
      Latitude = col_double(),
      Longitude = col_double()
    )
  )
}

#' Reads aru2point.csv
#'
#' This provides a manifest of ARU filenames and a mapping of filename to point
#' (where the ARU was at the time). Since the filename includes the time that
#' recording began, this is also is the authoritative manifest of ARU "visits."
#'
#' @return return tibble of (filename, point) for rows where point is defined.
#'
#' @export
readAru2point <- function() {
  read_csv(
    aru2pointPath,
    col_types = cols(
      filename = col_character(),
      point = col_integer(),
    )
  ) %>%
    filter(point > 0)
}

#' Filters to a time-of-day interval
#'
#' @param t tibble with a Date_Time column that informs the filtering.
#' @param beginTime Optional `lubridate::duration` to filter out rows with
#'   Date_Time less than this duration from their most recent midnight.
#' @param endTime Optional `lubridate::duration` to filter out rows with
#'   Date_Time greater than or equal to this duration from their most recent
#'   midnight.
#'
#' @return the filtered version of the input. This is intended to be used with
#'   the pipe operator.
filterTimeOfDay <- function(t, beginTime = NA, endTime = NA) {
  if (!is.na(beginTime)) {
    t <- t %>% filter((Date_Time - floor_date(Date_Time, "day")) >= beginTime)
  }
  if (!is.na(endTime)) {
    t <- t %>% filter((Date_Time - floor_date(Date_Time, "day")) < endTime)
  }
  t
}

#' Creates and populates a multidimensional array from (indices, values)
#'
#' This is similar to as.matrix(sparseMatrix), except that indices not included
#' in the entryTable param get a value of NA instead of 0.
#'
#' @param entryTable data.frame where the last column has the values to place
#'   into the returned array and where the rest of the columns, in order, are
#'   the indices where they will be placed.
#' @param indexTable data.frame with "index" columns corresponding to the names
#'   of all but the last column of entryTable. The columnwise maximums become
#'   the dimensions of the returned array.
#'
#' @return return array with the values-at-indices specified in entryTable and
#'   NA at indices where no value was specified.
sparseToDense <- function(entryTable, indexTable) {
  valueCol <- ncol(entryTable)
  indexNames <- names(entryTable)[-valueCol]
  dims <- sapply(indexTable[indexNames], max)
  
  indices <- entryTable[indexNames]
  values <- entryTable[, valueCol]
  
  dense <- array(NA, dims)
  dense[as.matrix(indices)] <- unlist(values)
  dense
}
