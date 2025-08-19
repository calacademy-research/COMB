# Code to iterate 7 different "knockout" model structures over a range of species

library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)

source(here("comb_functions.R"))
source(here("COMB_minimal/models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated# run "knockout" for each species on the list

set.seed(123)

# MCMC settings
na <- 100
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6

# parameters --------------------------------------------------------------
#speciesCode <- sort(c("BBWO", "WHWO", "NOFL", "HAWO", "PIWO", "WISA", "RBSA")) # must match prefiltering of dataML_model.csv
#speciesCode <- sort(c("AMRO", "WETA", "BHGR", "YRWA", "HETH", "WEWP", "SPTO"))
speciesCode <- sort(c("CONI", "BBWO", "AMRO", "WEWP", "RBNU", "OSFL", "MOUQ", "CORA", "CAHU", "HETH", "HEWA", "RBSA", "PUFI", "CAFI"))
speciesCode <- sort(c("SOGR", "ANHU"))
speciesCode <- sort(c("BBWO", "CONI"))

speciesCode <- sort(c("AMRO", "BBWO", "CONI", "DUFL"))
#year <- c(2020, 2021)
year <- 2021
threshold <- c(-2,-1,0,1,2)
aruVisitLimit <- 60 # only consider this many ARU visits per site (ordered). 32=4 full mornings


# read in/prepare field data ----------------------------------------------

# SITE covars
site_covars <- read_csv("COMB_minimal/models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

# JAGS structuring --------------------------------------------------------

for(th in threshold) {

for (sp in speciesCode) {
  begin_hour = 4
  end_hour = 10
  data <- readCombined(
    species = sp,
    years = c(year),
    beginTime = begin_hour,
    endTime = end_hour,
 #   visitLimit = aruVisitLimit, # total number of files
    visitAggregation = "file",
    daysLimit = 5,
    samplesperdayLimit = 12, # determines how many samples per day to sample
    thresholdOptions = list(value = th,
                            is.quantile = F),
    squeeze = T,
    logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
   # scale_datetime = T
  )

  covs <- site_covars %>% filter(Point %in% data$indices$point$Point) %>% arrange(Point)
  
  Veg <- left_join(as.data.frame(data$indices$point), covs, by = "Point")
  
  # Clean data object and add covariates
  jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',
                              "aru_DateTime","aru_YDay", "aru_Time","scores_datetime"))
  
  jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
  jagsData$burn  <- Veg$mean_RAVGcbi_20182019_4ha
  jagsData$Xburn <- seq(0, 2.5, length.out = 100)
  
  assign(paste("data", th, sp, year, sep="_"), jagsData)
}

datasets <- setNames(
  lapply(speciesCode, function(sp) get(paste("data", th, sp, year, sep="_"))),
  speciesCode
)


# monitored elements ------------------------------------------------------
# FULL (HUMAN + ARU + SCORES)
monitored_HAS <- c(
  "mean_psi",
  "beta0",
  "beta1",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "psi",
  "psi.pred.burn",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "TvegObs",
  "TvegSim",
  "Dobs",
  "Dsim"
)

#kos <- c("model_A", "model_AS", "model_H", "model_HA", "model_HAS", "model_HS", "model_S")
kos <- c("model_A", "model_AS", "model_HA", "model_HAS") # only ones we need for thresh expt
#kos <- "model_HAS"
#modDir <- "COMB_minimal/models/jags_files/"
modDir <- "COMB_minimal/models/jags_files/no_covars/"

#modDir <- "models/jags_files/"
#model <- paste0(modDir,"model_HAS_multiyear_nocovar.txt")
#sp <- "AMRO"
#data <- data_AMRO

# Inits lists -------------------------------------------------------------

results <- list()  # optional: store all results here

for (sp in names(datasets)) {
  data <- datasets[[sp]]
  
  # Create inits components specific to this data
  zst <- rep(1, data$nsites)
  psit <- runif(data$nsamples)
  
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > 0.0] <- 1
  gst[data$score <= 0.0] <- 2
  
  inits <- function() {
    list(
      mu = c(-2, 0),
      sigma = c(1, 1),
      z = zst,
      g = gst,
      beta0 = 0,
      beta1 = 0,
      reg_parm = 1
    )
  }

  for (k in seq_along(kos)) {
    model <- switch(kos[k],
                    "model_A"   = paste0(modDir, "model_A.txt"),
                    "model_AS"  = paste0(modDir, "model_AS.txt"),
                    "model_H"   = paste0(modDir, "model_H.txt"),
                    "model_HA"  = paste0(modDir, "model_HA.txt"),
                    "model_HAS" = paste0(modDir, "model_HAS.txt"),
                    "model_HS"  = paste0(modDir, "model_HS.txt"),
                    paste0(modDir, "model_S.txt"))  # default
    
    jagsResult <- jags(
      data = data,
      inits = inits,
      parameters.to.save = monitored_HAS,
      model.file = model,
      n.adapt = na,
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
    
    # for each result, give it a name and put it in the results list
    result_name <- paste("jagsResult", th, sp, paste(year, collapse = "-"),
                         str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
    
    assign(result_name, jagsResult)
    results[[result_name]] <- jagsResult
  }
#}
} # end of species dataset loop
rm(jagsResult)

jags_objects <- ls(pattern = "^jagsResult")

  # Loop through the list and save each object as an RData file
  for (obj in jags_objects) {
    save(list = obj, file = paste0("COMB_minimal/results/jagsResults/jagsResults_no_covars/", obj, ".RData"))
  }

}


# write summary table -----------------------------------------------------

params_of_interest <- c(
  "mean_psi", "p_11", "p_aru11", "p_aru01", 
  "mu[1]", "mu[2]", "sigma[1]", "sigma[2]", 
  "NOcc", "PropOcc", "deviance"
)

# Find all objects in the environment that start with "jagsResult_"
jags_objs <- mget(ls(pattern = "^jagsResult_"), inherits = TRUE)

# Function to process each model object
extract_summary <- function(obj, obj_name) {
  
  # Parse name: jagsResult_[threshold]_[species]_[year]_[model]_Sysdate
  parts <- unlist(strsplit(obj_name, "_", fixed = TRUE))
  
  tibble::as_tibble(obj$summary, rownames = "parameter") %>%
    filter(parameter %in% params_of_interest) %>%
    mutate(
      threshold = as.numeric(parts[2]),
      species   = parts[3],
      year      = as.integer(parts[4]),
      model     = paste(parts[6:(length(parts) - 1)], collapse = "_"),
      sysdate   = parts[length(parts)]
    )
}

# Apply to all objects and combine
df_all <- imap_dfr(jags_objs, extract_summary)

# Check result
print(df_all)

write.csv(df_all, file = paste0("COMB_minimal/results/jagsResults/jagsResults_no_covars/threshold_expt_3spp", Sys.Date(), ".csv"))

