#Spatial-temporal leave one out cross validation

#######################
# Backwards variable selection using sloo.loglik on full model
#######################

# Arguments:
# - full model
# - training data sets from splitting function
modSelect <- function(full.mod, data, training, steps){
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else extractAIC(x, k = 0)[2L]
  }
  object <- full.mod
  Terms <- terms(object)
  # object$call$formula <- object$formula <- Terms
  fdrop <- numeric()
  fadd <- attr(Terms, "factors")
  models <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- fitmod <- object
  # nfit <- fit$call
  bAIC <- extractAIC(fit)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  # bloglik <- stloo.loglik(model = fit, training = training, nfit = nfit)
  # bRMSE <- stloo.rmse(model = fit, training = training, data = data)
  bfit <- fullmodfit <- stCV.rmse(model = fit, training = training, data = data)
  bRMSE <- bfit$rmse
  bR2train <- bfit$meanR2train
  bR2test <- bfit$meanR2test
  
  
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'step' cannot proceed")
  if (bAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'step' cannot proceed")
  nm <- 1
  # models[[nm]] <- list(fit = fit, deviance = mydeviance(fit), df.resid = n - 
  #                        edf, change = "", AIC = bAIC, loglik = bloglik)
  models[[nm]] <- list(fit = fitmod, deviance = mydeviance(fit), change = "", AIC = bAIC, rmse = bRMSE, 
                       meanR2test = bR2test, meanR2train = bR2train, drops = 0, 
                       fullrmse = fullmodfit$rmse, fullmeanR2test = fullmodfit$meanR2test, fullmeanR2train = fullmodfit$meanR2train)
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    # loglik <- bloglik
    rmse <- bRMSE
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (length(scope$drop)) {
      aod <- drop1stloo(object = fitmod, scope = scope$drop, training = training, data = data)
      rn <- row.names(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df)) 
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0) 
        break
      nc <- match(c("rmse", "loglik", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc], decreasing = FALSE)
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    
    fitmod <- update(fitmod, paste("~ .", change), control = list(maxIter = 100, opt = "optim"))
    # fitmod <- eval.parent(fit)

    Terms <- terms(fitmod)
    bAIC <- extractAIC(fitmod)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    # bloglik <- stloo.loglik(model = fit, training = training, nfit = fit)
    # bRMSE <- stloo.rmse(model = fitmod, training = training, data = data)
    bfit <- stCV.rmse(model = fitmod, training = training, data = data)
    bRMSE <- bfit$rmse
    bR2test <- bfit$meanR2test
    bR2train <- bfit$meanR2train
    # if (bloglik <= loglik - 1e-05) 
    #   break    

    if (bRMSE >= rmse + 1e-05) 
      break  
    
    nm <- nm + 1
    # models[[nm]] <- list(fit = fit, deviance = mydeviance(fit), df.resid = n - 
    #                        edf, change = change, AIC = bAIC, loglik = bloglik, drops = nm)
    models[[nm]] <- list(fit = fitmod, deviance = mydeviance(fitmod), 
                         change = change, AIC = bAIC, rmse = bRMSE, meanR2test = bR2test,
                         meanR2train = bR2train, drops = (nm - 1),
                         fullrmse = fullmodfit$rmse, fullmeanR2test = fullmodfit$meanR2test, fullmeanR2train = fullmodfit$meanR2train)
  }
 
return(models[[nm]])  
}


drop1stloo <- function (object, scope, training, data) 
{
  tl <- attr(terms(object), "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object)
  }else {
    if (!is.character(scope)) {
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    }
    if (!all(match(scope, tl, 0L) > 0L)) {
      stop("scope is not a subset of term labels")
    }
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                             scope), c("df", "rmse")))
  n <- length(object$residuals)
  if(is.null(object$df.residaul)){
    object$df.residual <- object$dims$N - object$dims$p
  }
  edf <- n - object$df.residual
  # ans[1, ] <- c(edf, stloo.loglik(object, training, object))
  # ans[1, ] <- c(edf, stloo.rmse(object, training, data))
  ans[1, ] <- c(edf, stCV.rmse(object, training, data)$rmse)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq_len(ns)) {
    tt <- scope[i]
    
    nfit <- update(object, as.formula(paste("~ . -", tt)), 
                   control = list(maxIter = 100, opt = "optim"))
    
    # ans[i + 1, ] <- c(edf + 1, stloo.loglik(object, training, nfit))
    # ans[i + 1, ] <- c(edf + 1, stloo.rmse(nfit, training, data))
    ans[i + 1, ] <- c(edf + 1, stCV.rmse(nfit, training, data)$rmse)
  }
  dfs <- ans[1L, 1L] - ans[, 1L]
  dfs[1L] <- NA
  # aod <- data.frame(Df = dfs, loglik = ans[, 2])
  aod <- data.frame(Df = dfs, rmse = ans[, 2])
  aod
}

#########################################################################
# Splitting the data from its distance matrix and a threshold distance ##
#########################################################################

# Arguments:
# - data: the 'data.frame' used
# - lon.lab, lat.lab: character strings giving the labels of longitude and latitude columns in the data
# - dist.split: the threshold distance

splittingST <- function(data, lon.lab, lat.lab, yr.lab, dist.split, time.split){
  # Computing the distance matrix from the data:
  # dist.matrix <- as.matrix(x=dist(x=data[,c(lon.lab,lat.lab)], method = "euclidean",diag = T, upper = T))
  dist.matrix <- GeoDistanceInMetresMatrix(data.frame(lat = data[, lat.lab], lon = data[, lon.lab]))
  time.matrix <- as.matrix(x=dist(x=data[, yr.lab], method = "euclidean",diag = T, upper = T))
  
  # Initializing the 'training' object
  training <- list()
  # Creating the training sets
  for(i in 1:nrow(dist.matrix)){
    # Keeping only the observations far enough of the i-st observation by using the threshold distance
    num.cell <- as.vector(which(dist.matrix[i,] > dist.split & time.matrix[i, ] > time.split))
    training[[i]] <- num.cell
  }
  return(training)
}

# split data via ENMeval package get.block method
# after first separating train and test sets by time
splittingCV <- function(data, lon.lab, lat.lab, yr.lab, time.split){
  data$index <- 1:nrow(data)
  data$timelevel <- Hmisc::cut2(data[, yr.lab], g = time.split)
  
  allpts <- data[ , c("lon", "lat", "SiteID", "index")]
  
  allptslat <- Hmisc::cut2(allpts[, lat.lab], g = 2)
  
  # modified from ENMeval::get.block
  grpA <- allpts[which(allptslat == levels(allptslat)[1]), ]
  grpB <- allpts[which(allptslat == levels(allptslat)[2]), ]
  grpAlon <- Hmisc::cut2(grpA[, lon.lab], g = 2)
  grpBlon <- Hmisc::cut2(grpB[, lon.lab], g = 2)
  
  grp1 <- grpA[which(grpAlon == levels(grpAlon)[1]), ]
  grp2 <- grpA[which(grpAlon == levels(grpAlon)[2]), ]
  grp3 <- grpB[which(grpBlon == levels(grpBlon)[1]), ]
  grp4 <- grpB[which(grpBlon == levels(grpBlon)[2]), ]
  grp1$geolevel <- "SW"
  grp2$geolevel <- "SE"
  grp3$geolevel <- "NW"
  grp4$geolevel <- "NE"
  grp <- rbind(grp1, grp2, grp3, grp4)
  
  grp <- merge(grp, data[, c("SiteID", "index", "timelevel")], by = c("SiteID", "index"))
  
  # selecting for training data: 
  # all regions for years not removed to test
  # 3 regions for years removed to test
  # amounts roughly to 12-fold CV
  training <- list()
  for(reg in unique(grp$geolevel)){
    for (yr in unique(grp$timelevel)){
      trainrows <- grp$index[-which(grp$geolevel == reg & grp$timelevel == yr)]
      num.cell <- as.vector(sort(trainrows))
      if(length(num.cell) >= (nrow(data) - 2)) next  #does not create dataset if <2 test pts
      training[[length(training) + 1]] <- num.cell
    }
  }
  return(training)
}
######################################################################
# Computing Spatial Leave One Out (SLOO) logLikelihood from a model ##
######################################################################

# Arguments:
# - model: an object of class 'lm' or 'glm' giving the model considered.
# - training: the 'list' return by the 'splitting' function (see splitting.r).

stloo.rmse <- function(model, training, data){
  modeldat <- as.data.frame(data)
  predmod <<- model
  predform <<- formula(predmod)
  # Creating the response variable name
  y <- as.character(x=formula(x=predmod)[2])
  # Initializing the 'logLik' object
  sqerr <- vector(mode="numeric",length=nrow(modeldat))
  # Checking that the minimal number of observations in the training sets
  # is not lower than the number of parameters in the model considered.
  if(min(as.numeric(x=summary(object=training)[,1]))<(length(x=predmod$coef)+1)){
    # Return NA if the minimal number of observations in the training sets
    # is lower than the number of parameters of the model considered
    print(x="Warning: too high threshold distance in 'training', NA is return")
    return(value=NA)
  }
  else {
    # Calculating the SLOO logLikelihoods for each observation i
    for(i in 1:nrow(modeldat)){ 
      # Extracting the i-st training set:
      data.split <- as.data.frame(modeldat[training[[i]],])
      # Calculating the model parameters from the i-st training set:
      # m <- glm(formula=formula(nfit),data=data.split,family= model$family)
      m <- gls(model = predform, 
               data = data.split, correlation = corExp(form = ~ lat + lon | YearFact), method = "ML")
      # Predicting the i-st observation from the i-st training set:
      m.pred <- predict(object=m, newdata=modeldat[i,])
      # Calculating the probability of the i-st observed value according to the predicted one by the i-st training set:
      sqerr[i] <- (modeldat[i,y] - m.pred)^2
      
    }
    # Calculating the overall SLOO rmse:
    rmse <- sqrt(mean(sqerr))
    return(value = rmse)
  }
}

# LOO too computationally expensive
# this function replaces LOO with approx 12-fold cross validation
# time split observations into 3 equal groups, each group split into 4 regions
stCV.rmse <- function(model, training, data){
  modeldat <- as.data.frame(data)
  predmod <<- model
  predform <<- formula(predmod)
  # Creating the response variable name
  y <- as.character(x=formula(x=predmod))[2]
  nullform <<- formula(paste(y, "~ 1", sep = " "))

  # Initializing the 'logLik' object
  sqerr <- rep(NA,length(training))
  r2 <- rep(NA,length(training))
  r2train <- rep(NA, length(training))
  r2test <- rep(NA, length(training))

  
  # Checking that the minimal number of observations in the training sets
  # is not lower than the number of parameters in the model considered.
  if(min(as.numeric(x=summary(object=training)[,1]))<(length(x=predmod$coef)+1)){
    # Return NA if the minimal number of observations in the training sets
    # is lower than the number of parameters of the model considered
    print(x="Warning: too high threshold distance in 'training', NA is return")
    return(value=NA)
  }else {
    # Calculating the SLOO logLikelihoods for each observation i
    for(i in 1:length(training)){ 
      # Extracting the i-st training set:
      data.train <- as.data.frame(modeldat[training[[i]], ])
      data.test <- as.data.frame(modeldat[-training[[i]], ])
      # Calculating the model parameters from the i-st training set:
      # m <- glm(formula=formula(nfit),data=data.split,family= model$family)
      m <- try(gls(model = predform, 
                   data = data.train, 
                   correlation = corExp(form = ~ lat + lon | YearFact),
                   method = "ML",
                   control = list(maxIter = 100, opt = "optim"))
      )
      null <- try(gls(model = nullform,
                      data = data.train,
                      correlation = corExp(form = ~ lat + lon | YearFact),
                      method = "ML",
                      control = list(maxIter = 100, opt = "optim"))
      )
      if (class(m) == "gls"){
        # Predicting the i-st observation from the i-st training set:
        m.pred <- predict.gls.tyson(object = m, newdata = data.test)
        null.pred <- predict.gls.tyson(object = null, newdata = data.test)
        # insample and out of sample nagelkerke r2
        m.loglik <- sum(dnorm(x = as.numeric(data.test[,y]), mean=m.pred, sd=sqrt(sum(residuals(m)^2)/nrow(data.train)),log=TRUE))
        null.loglik <- sum(dnorm(x = as.numeric(data.test[,y]), mean=null.pred, sd=sqrt(sum(residuals(null)^2)/nrow(data.train)),log=TRUE))
        r2train[i] <- nagR2(m$dims$N, as.numeric(logLik(m)), as.numeric(logLik(null)))
        r2test[i] <- nagR2(length(m.pred), m.loglik, null.loglik)
        # Calculating the probability of the i-st observed value according to the predicted one by the i-st training set:
        sqerr[i] <- mean((data.test[,y] - m.pred)^2)
        # r2[i] <- 1 - sum((data.test[,y] - m.pred)^2) / sum((data.test[,y] - mean(data.test[,y]))^2)
      }else{
        next
      }
    }
    # Calculating the overall SLOO rmse:
    rmse <- sqrt(mean(sqerr))
    
    # meanR2 <- mean(r2)
    meanR2train <- mean(r2train)
    meanR2test <- mean(r2test)
    return(data.frame(rmse = rmse, meanR2train = meanR2train, meanR2test = meanR2test))
  }
}

nagR2 <- function(nobs, ll.full, ll.null){
  r2 <- 1 - exp((-2/nobs) * (ll.full - ll.null))
  r2max <- 1 - exp((2/nobs) * ll.null)
  r2corr <- r2/r2max
  return(r2corr)
}

nagR2v2 <- function(mod.full, mod.null){
  ll.full <- as.numeric(logLik(mod.full))
  ll.null <- as.numeric(logLik(mod.null))
  nobs <- mod.full$dims$N
  if (nobs == mod.null$dims$N){
  r2 <- 1 - exp((-2/nobs) * (ll.full - ll.null))
  r2max <- 1 - exp((2/nobs) * ll.null)
  r2corr <- r2/r2max
  return(r2corr)
  }else{
    print("Error!: Dimensions don't match")
  }
}
######################################################################
# Computing Spatial Leave One Out (SLOO) logLikelihood from a model ##
######################################################################

# Arguments:
# - model: an object of class 'lm' or 'glm' giving the model considered.
# - training: the 'list' return by the 'splitting' function (see splitting.r).

stloo.loglik <- function(model, training, nfit){
  modeldat <- as.data.frame(model$data)
  training <- training
  # Creating the response variable name
  y <- as.character(x=formula(x=model)[2])
  # Initializing the 'logLik' object
  logLik <- vector(mode="numeric",length=nrow(x=model$data))
  # Checking that the minimal number of observations in the training sets
  # is not lower than the number of parameters in the model considered.
  if(min(as.numeric(x=summary(object=training)[,1]))<(length(x=model$coef)+1)){
    # Return NA if the minimal number of observations in the training sets
    # is lower than the number of parameters of the model considered
    print(x="Warning: too high threshold distance in 'training', NA is return")
    return(value=NA)
  }
  else {
    # Calculating the SLOO logLikelihoods for each observation i
    for(i in 1:nrow(x=model$data)){ 
      # Extracting the i-st training set:
      data.split <- as.data.frame(modeldat[training[[i]],])
      # Calculating the model parameters from the i-st training set:
      m <- glm(formula=formula(nfit),data=data.split,family= model$family)
      # Predicting the i-st observation from the i-st training set:
      m.pred <- predict(object=m,newdata=modeldat[i,],type="response")
      # Calculating the probability of the i-st observed value according to the predicted one by the i-st training set:
      if(model$family[1]=="gaussian") logLik[i] <- dnorm(x = as.numeric(modeldat[i,y]), mean=m.pred, sd=sqrt(sum(residuals(m)^2)/nrow(data.split)),log=TRUE)
      else if(model$family[1]=="binomial"){
        if(length(unique(model$y))!=2 | !is.element(0,unique(model$y)) | !is.element(1,unique(model$y))){
          print("The calculation of the Spatial Leave One Out logLikelihood") 
          print("could not be done by this function for non-binary response variable.")
          print("Please reconsider the response variable of the model in 0/1 coding scheme.")
          break
        }  
        else logLik[i] <- dbinom(x=modeldat[i,y],size=1,prob=m.pred,log=T)
      }  
      else if(model$family[1]=="poisson")  logLik[i] <- dpois(x=modeldat[i,y],lambda=m.pred,log=T)
      else{ 
        print("This family is not supported by this function")
        break
      }
    }
    # Calculating the overall SLOO logLikelihood:
    Sum.logLik <- sum(logLik)
    return(value=Sum.logLik)
  }
}


# from http://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.

  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.

  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }

  n.geopoints <- nrow(df.geopoints)

  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints

  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})

  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")

  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name

  return(mat.distances)
}

# EDIT predict.gls so that factor levels work 
# see http://stackoverflow.com/questions/22229788/why-does-predict-on-a-subset-of-the-original-data-fail-using-gls-when-lm-does

predict.gls.tyson <- function (object, newdata, na.action = na.fail, ...) 
{
  if (missing(newdata)) {
    return(fitted(object))
  }
  form <- nlme::getCovariateFormula(object)
  mfArgs <- list(formula = form, data = newdata, na.action = na.action)
  # mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call(model.frame, mfArgs)
  contr <- object$contrasts
  for (i in names(dataMod)) {
    if (inherits(dataMod[, i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMod[, i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s", 
                              "levels %s not allowed for %s"), paste(levs[wch], 
                                                                     collapse = ",")), domain = NA)
      }
      attr(dataMod[, i], "contrasts") <- contr[[i]][levs, 
                                                    , drop = FALSE]
    }
  }
  N <- nrow(dataMod)
  if (length(all.vars(form)) > 0) {
    X <- model.matrix(form, dataMod)
  }
  else {
    X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
  }
  cf <- coef(object)
  val <- c(X[, names(cf), drop = FALSE] %*% cf)
  lab <- "Predicted values"
  if (!is.null(aux <- attr(object, "units")$y)) {
    lab <- paste(lab, aux)
  }
  structure(val, label = lab)
}