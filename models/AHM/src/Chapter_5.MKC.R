# Chapter 5 in AHM2 introdues a dynamic community model or DCM (multi-year, multi-species).

rm(list = ls())

library(AHMbook)


# 5.2: A general simulation function for the DCM --------------------------



# play with how you can manipulate the data generation process...
# Explicit defaults (need to load AHMbook)
# str(data <- simDCM(nspec = 50, nsites = 100, nsurveys = 3, nyears = 10,
#                    mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
#                    range.mean.phi = c(0.8, 0.8), sig.lphi = 1, mu.beta.lphi = 0,
#                    sig.beta.lphi = 0, range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
#                    mu.beta.lgamma = 0, sig.beta.lgamma = 0, range.mean.p = c(0.5, 0.5),
#                    sig.lp = 1, mu.beta.lp = 0, sig.beta.lp = 0, range.beta1.survey = c(0, 0),
#                    range.beta2.survey = c(0, 0), trend.sd.site = c(0, 0),
#                    trend.sd.survey = c(0, 0), show.plot = TRUE) )
# str(data <- simDCM(nspec = 200)) # More species (looks great) str(data <- simDCM(nspec = 1)) # A single species (works !) str(data <- simDCM(nsites = 267)) # More sites
# str(data <- simDCM(nsites = 1)) # A single site
# str(data <- simDCM(nsurveys = 10)) # More visits
# str(data <- simDCM(nyears = 25))
# str(data <- simDCM(nyears = 2))
# try(data <- simDCM(nyears = 1))
# # More years
# # Just two years
# # A single year ... this crashes
# # No species heterogeneity in parameters of initial occupancy
# str(data <- simDCM(sig.lpsi1 = 0, sig.beta.lpsi1 = 0))
# # No species heterogeneity in parameters of persistence
# str(data <- simDCM(sig.lphi = 0, sig.beta.lphi = 0))
# # No species heterogeneity in parameters of colonization
# str(data <- simDCM(sig.lgamma = 0, sig.beta.lgamma = 0))
# # No species heterogeneity in parameters of detection
# str(data <- simDCM(sig.lp = 0, sig.beta.lp = 0))
# # No annual variation in rates phi, gamma and p
# str(data <- simDCM(range.mean.phi = c(0.8, 0.8), range.mean.gamma = c(0.3, 0.3),
#                    range.mean.p = c(0.6, 0.6)))


# 5.4: Fun with Multidimensional Arrays -----------------------------------

# We illustrate the following steps: (1) We take a 4D array and shoot “holes” by turning 25% of the data into missing values, (2) we vectorize the data, (3) we get rid of the NAs, resulting in a data format as if it came to you in an Excel spreadsheet, (4) and then we re-assemble the spreadsheet-format data back into a 4D array.

# Get a 4D array

set.seed(123)
tmp <- simDCM(show.plot = F)
BigArray <- tmp$y # Grab the 4D detection/nondetection array str(BigArray) # 100 sites x 3 visits x 10 years x 50 species

# (1) Shoot holes, i.e., randomly turn 25% of the values into NAs
length(BigArray) # 150k values
out <- sample(1:length(BigArray), length(BigArray) / 4) # Sample of 25%
BigArray[out] <- NA # Turn them into NAs
sum(is.na(BigArray)) # we now have 37500 NAs

# (2) Turn data into vector and create indexing factors for the array dimensions
df <- data.frame(
  y = c(BigArray),
  site = c(slice.index(BigArray, 1)), visit = c(slice.index(BigArray, 2)), year = c(slice.index(BigArray, 3)), species = c(slice.index(BigArray, 4))
)

summary(df) # Check the 'max.' for each column
str(df)

# (3) Toss out the rows with missing response
df <- df[!is.na(df$y), ] # Select rows with non-missing y only
sum(is.na(df$y)) # Convince yourself NAs are gone
head(df) # Look at first 6 rows in data frame

# (4) Format data in a spreadsheet format into a 4D array # Determine required dimensions of 4D array
nsite <- length(unique(df$site))
nvisit <- length(unique(df$visit))
nyear <- length(unique(df$year))
nspec <- length(unique(df$species))
# Prepare array and pre-fill array with NAs
# Number of sites
# Number of surveys or visits
# Number of years
# Number of species
BigArray2 <- array(NA, dim = c(nsite, nvisit, nyear, nspec))

# Fill array with the detection/nondetection data
# Loop over all rows in the spreadsheet data and fill them in # at the right place in the 4D array
for (i in 1:nrow(df)) {
  BigArray2[df$site[i], df$visit[i], df$year[i], df$species[i]] <- df$y[i]
}
# Do quick checks ... look good
sum(df$y)
sum(BigArray2, na.rm = TRUE) # quick sum check
length(which(is.na(BigArray2))) # Same 37500 as before
all.equal(BigArray, BigArray2, check.attributes = FALSE) # BigArray has names

#
# Get a data set with detections only, no nondetections (similar to our data)
str(BigArray <- tmp$y) # Grab the 4D detection/nondetection array again
df <- data.frame(
  y = c(BigArray),
  site = c(slice.index(BigArray, 1)),
  visit = c(slice.index(BigArray, 2)),
  year = c(slice.index(BigArray, 3)),
  species = c(slice.index(BigArray, 4))
) # Turn into a spreadsheet format again
df <- df[df$y == 1, ] # Toss out all nondetection data

# Prepare array by pre-filling it with zeroes instead of NAs
BigArray3 <- array(0, dim = c(100, 3, 10, 50)) # known array dims !
# Fill array with the detection data
for (i in 1:nrow(df)) {
  BigArray3[df$site[i], df$visit[i], df$year[i], df$species[i]] <- df$y[i]
}
sum(BigArray3) - nrow(df) # quick check they're identical

# Reformat array so sites come last
dim(BigArray3) # 100   3  10  50
dim(BA.v2 <- aperm(BigArray3, c(2, 3, 4, 1))) # 3  10  50 100

head(df)
