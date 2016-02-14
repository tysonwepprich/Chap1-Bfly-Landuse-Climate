# Chapter 1 functions

#Outlier detection, problem with weird sightings changing GAMs
OutlierDetect <- function(counts){
  counts$DaysOut <- NA
  for (i in 1:nrow(counts)){
    allelse <- counts[-i]
    temp <- allelse[Year == counts[i, Year]][Region == counts[i, Region]]
    if (counts[i, Total] == 0){
      counts$DaysOut[i] <- 0
    }else{
      if (length(temp$Ordinal[temp$Total > 0]) == 0){
        counts$DaysOut[i] <- NA
      }else{
        counts$DaysOut[i] <- min(abs(counts[i, Ordinal] - temp$Ordinal[temp$Total > 0]))
      }
    }
  }
  return(counts)
}


#Then Dennis gets a little confusing
#Equation 4 needs the scaled phenology for each year and region from the GAM predictions
#Then fits Poisson GLM to counts with site and year as (non-interacting) factors
#Scaled phenology is used as an offset in the GLM
ScaleSumTo1 <- function(x){x/sum(x)}

ScaleDennis <- function(x){exp(x)/sum(exp(x))}

ScaledPhenology <- function(counts, yr, r){
  temp <- counts[Year == yr][Region == r]
  knts <- 10 # try this until it breaks down
  #   if (dim(temp)[1] < 10){
  #     knts <- dim(temp)[1]
  #   }else{
  #     knts <- 10
  #   }
  
  if(length(unique(temp$SiteID)) == 1){
    species.model <- gam(Total ~ s(Ordinal, bs = "cr", k = knts), 
                         family = negbin(theta = c(1,10), link = "log"), data = temp, 
                         method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    
    start.surv <- min(temp$Ordinal)
    end.surv <- max(temp$Ordinal)
    #     start.surv <- 91
    #     end.surv <- 305
    
    GAM.pred <- predict.gam(species.model, data.frame(Ordinal = c(start.surv:end.surv)), type = "response")
    output <- data.frame(Year = yr, Region = r, Ordinal = c(start.surv:end.surv), Gamma = ScaleSumTo1(GAM.pred))
  }else{
    
    species.model <- gam(Total ~ SiteID + s(Ordinal, bs = "cr", k = knts), 
                         family = negbin(theta = c(1,10), link = "log"), data = temp, 
                         method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    
    start.surv <- min(temp$Ordinal)
    end.surv <- max(temp$Ordinal)
    #       start.surv <- 91
    #       end.surv <- 305
    maxSite <- temp$SiteID[temp$YearTotal == max(temp$YearTotal)][1] #doesn't really matter which site
    
    GAM.pred <- predict.gam(species.model, data.frame(Ordinal = c(start.surv:end.surv), SiteID = maxSite), type = "response")
    output <- data.frame(Year = yr, Region = r, Ordinal = c(start.surv:end.surv), Gamma = ScaleSumTo1(GAM.pred))
  }
  return(output)
}


ScaledPhenologyNB <- function(counts, yr){
  temp <- counts[Year == yr]
  # add zero counts a month outside of monitoring season
  zeroCounts <- expand.grid(unique(temp$SiteID), c(59,60, 334,335))
  names(zeroCounts) <- c("SiteID", "Ordinal")
  zeroCounts$Total <- 0
  
  geo <- unique(temp[, c("SiteID", "lat", "lon"), with = FALSE])
  
  zeroCounts <- merge(zeroCounts, geo, by = "SiteID")
  temp <- rbind(temp, zeroCounts, fill = TRUE)
  #   temp[, AvgSitePop := mean(log(Total + 1/30)), by = SiteID]
  knts <- c(20, min(c(5, length(unique(temp$SiteID)) - 1))) # try this until it breaks down
  
  if (max(temp$lat) - min(temp$lat) < 1){
    mod <- gam(Total ~ SiteID + s(Ordinal, k = knts[1]),
               family = nb(theta = NULL, link = "log"), data = temp, 
               method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
  }else{
    
    mod <- gam(Total ~ SiteID + te(Ordinal, lat, k = knts),
               family = nb(theta = NULL, link = "log"), data = temp, 
               method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
    #       
    #           mod <- gam(Total ~ SiteID + te(Ordinal, lat, k = knts),
    #                  family = negbin(theta = prelim$family$getTheta(TRUE), link = "log"), data = temp, 
    #                  method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))
  }
  
  # get actual start and end, without padded zeros
  start.surv <- min(sort(unique(temp$Ordinal))[3:(length(unique(temp$Ordinal)) - 2)])
  end.surv <- max(sort(unique(temp$Ordinal))[3:(length(unique(temp$Ordinal)) - 2)])
  
  maxSite <- temp$SiteID[temp$YearTotal == max(temp$YearTotal)][1] #doesn't really matter which site
  
  newData <- expand.grid(unique(temp$SiteID), c(start.surv:end.surv))
  names(newData) <- c("SiteID", "Ordinal")
  newData$Year  <- yr
  newData <- merge(newData, geo, by = "SiteID")
  
  newData$GAM.pred <- predict.gam(mod, newData, type = "response")
  newData <- data.table(newData)
  newData[, Gamma := ScaleSumTo1(GAM.pred), by = "SiteID"]
  
  # plotting to check results
#   library(ggplot2)
  # a <- ggplot(data = newData, aes(x = Ordinal, y = GAM.pred, group = SiteID, color = SiteID)) + geom_line()
#   
  return(newData)
}


FillbyWeek <- function(counts){
  years <- sort(as.numeric(unique(counts$Year)))
  out <- data.frame()
  for (year in years){
    temp <- counts[Year == year]
    
    AllWeeks <- data.frame(expand.grid(c(unique(temp$SiteID)), c(1:max(temp$Week))))
    names(AllWeeks) <- c("SiteID", "Week")
    AllWeeks$Total <- NA
    
    if(year %% 4 == 0){     #Leap year changes ordinal dates
      AllWeeks$Ordinal <- 92 + (AllWeeks$Week - 1) * 7 + 3 #Pick middle day of week to predict missing count
    }else{
      AllWeeks$Ordinal <- 91 + (AllWeeks$Week - 1) * 7 + 3
    }
    AllWeeks$Year <- year
    out <- rbind(out, AllWeeks)
  }
  return(out)
}

MissingDays <- function(counts){
  years <- sort(as.numeric(unique(counts$Year)))
  out <- data.frame()
  for (year in years){
    temp <- counts[Year == year]
    
    AllWeeks <- data.frame(expand.grid(c(unique(temp$SiteID)), c(1:max(temp$Week))))
    names(AllWeeks) <- c("SiteID", "Week")
    AllWeeks$Total <- NA
    
    test <- merge(temp, AllWeeks, all.y = TRUE, by = c("SiteID", "Week"))
    missing_weeks <- test[is.na(Total.x) == TRUE, list(SiteID, Week)]
    if(year %% 4 == 0){     #Leap year changes ordinal dates
      missing_weeks$Ordinal <- 92 + (missing_weeks$Week - 1) * 7 + 3 #Pick middle day of week to predict missing count
    }else{
      missing_weeks$Ordinal <- 91 + (missing_weeks$Week - 1) * 7 + 3
    }
    missing_weeks$Year <- year
    out <- rbind(out, missing_weeks)
  }
  return(out)
}

#From UKBMS method, Equation 2 in Dennis, trapezoidal rule

TrapezoidIndex <- function(t, y){
  dat <- data.frame(cbind(t, y))
  dat <- arrange(dat, t)
  temp <- numeric()
  for (i in 2:length(dat$t)){
    temp[i-1] <- (dat$y[i] + dat$y[i-1]) * (dat$t[i] - dat$t[i-1]) / 2
  }
  return(sum(temp))
}