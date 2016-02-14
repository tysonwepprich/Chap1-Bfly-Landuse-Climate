# data prep for Bayesian models in jags
# source this file before running models
# produces multispecies matrix of counts, covariate array

# TODO: wrap all this into cleaner function
# TODO: add additional covariates for landuse at sites

library(plyr)
library(dplyr)
library(data.table)
library(reshape2)
library(rjags)
library(R2jags)


checklist <- read.csv("data/Checklist.csv", header = TRUE)
checklist <- checklist[, c("CheckListKey", "CommonName")]
names(checklist)[1] <- "Species"
checklist <- rbind(checklist, data.frame(Species = 4363, CommonName = "Azures"))
load("data/all.species.data.RData")
names(all.species.env)[c(1,2,5)] <- c("SiteID", "Year", "Species")

data <- readRDS(file = "data/AllSpeciesPops.rds")
data <- merge(data, checklist, by = "CommonName", all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE)
data <- merge(data, all.species.env, by = c("SiteID", "Year", "Species"))

JAGSdat <- data[, c(1, 2, 4, 5, 9, 15, 16), with = FALSE]
JAGSdat$TrpzInd <- round(JAGSdat$TrpzInd)
JAGSdat$TrpzNt1 <- round(JAGSdat$TrpzNt1)

# SpeciesList <- readRDS(file = "data/SpeciesList.rds")
# species <- as.character(SpeciesList$CommonName)
# species[4] <- "Azures"
# species <- species[1:68] 
# 
# # choose species
# spec_sel <- SpeciesList[c(10:20), "CommonName"]
# JAGSdat <- JAGSdat %>%
#   filter(CommonName %in% spec_sel)

dat <- JAGSdat %>%
  group_by(SiteID, CommonName) %>%
  mutate(LengthTS = max(Year) - min(Year), NumCounts = length(unique(Year))) %>%
  filter(NumCounts >= 5)

# what to do about missing years? 
# can estimate in time series approach, but think I need the missing covariates, too
dat$Year <- as.factor(as.character(dat$Year))
dat[, SiteID := formatC(SiteID, width = 3, format = "d", flag = "0")]

# all years x sites
allsiteyr <- expand.grid(sort(unique(dat$Year)), sort(unique(dat$SiteID)), sort(unique(dat$CommonName)))
names(allsiteyr) <- c("Year", "SiteID", "CommonName")
alldat <- merge(allsiteyr, dat, by = c("Year", "SiteID", "CommonName"), all.x = TRUE)
# assign zeros when surveyed but not counted, NA when not surveyed that year/site
for (i in 1:nrow(alldat)){
  temp <- alldat[i, ]
  if (is.na(temp$TrpzInd) == TRUE){
    site_yr <- alldat %>% filter(Year == temp$Year & SiteID == temp$SiteID)
    site_yr <- site_yr[complete.cases(site_yr), ]
    if (nrow(site_yr) > 0){
      alldat$TrpzInd[i] <- 0
      alldat$win.Comp.1[i] <- site_yr$win.Comp.1[1]
      alldat$win.Comp.2[i] <- site_yr$win.Comp.2[1]
    }
  }
}



# model covariates and predict missing values
env.dat <- unique(dat[, c("SiteID", "Year", "win.Comp.1", "win.Comp.2"), with = FALSE])

cov.mod.1 <- lm(win.Comp.1 ~ SiteID + Year, data = env.dat)
cov.mod.2 <- lm(win.Comp.2 ~ SiteID + Year, data = env.dat)

alldat[is.na(alldat$win.Comp.1), "win.Comp.1"] <- predict(cov.mod.1, newdata = alldat[is.na(alldat$win.Comp.1), ] )
alldat[is.na(alldat$win.Comp.2), "win.Comp.2"] <- predict(cov.mod.2, newdata = alldat[is.na(alldat$win.Comp.2), ] )

# cut out 1997, will use priors for this years starting population sizes
finaldat <- alldat[-which(alldat$Year == "1997"), 1:7]

nyears <- length(unique(finaldat$Year))
nsites <- length(unique(finaldat$SiteID))
nspecies <- length(unique(finaldat$CommonName))

# put all data into an array for JAGS (not necesary, but seems organized like Kery & Schaub)
count_array <- array(NA, dim = c(nsites, nyears, nspecies))
cov_array <- array(NA, dim = c(nsites, nyears, 2))
cov_molten <- melt(finaldat, id = c("SiteID", "Year", "CommonName"))

count_array <- acast(cov_molten[cov_molten$variable == "TrpzInd", ], SiteID ~ Year ~ CommonName, value = "value")
cov_array[,,1] <- acast(cov_molten[cov_molten$variable == "win.Comp.1", ], SiteID ~ Year ~ CommonName, value = "value")[,,1]
cov_array[,,2] <- acast(cov_molten[cov_molten$variable == "win.Comp.2", ], SiteID ~ Year ~ CommonName, value = "value")[,,1]

Nstartvals <- round(apply(count_array, c(1,3), mean, na.rm = TRUE))
