#Spatial-temporal leave one out cross validation

#######################
# Backwards variable selection using sloo.loglik on full model
#######################

# Arguments:
# - full model
# - training data sets from splitting function
modSelect <- function(full.mod, training, steps){
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev)) 
      dev
    else extractAIC(x, k = 0)[2L]
  }
  training <- training
  object <- full.mod
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  fdrop <- numeric()
  fadd <- attr(Terms, "factors")
  models <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  nfit <- fit$call
  bAIC <- extractAIC(fit)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  bloglik <- stloo.loglik(model = fit, training = training, nfit = nfit)
  if (is.na(bAIC)) 
    stop("AIC is not defined for this model, so 'step' cannot proceed")
  if (bAIC == -Inf) 
    stop("AIC is -infinity for this model, so 'step' cannot proceed")
  nm <- 1
  models[[nm]] <- list(fit = fit, deviance = mydeviance(fit), df.resid = n - 
                         edf, change = "", AIC = bAIC, loglik = bloglik)
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    loglik <- bloglik
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (length(scope$drop)) {
      aod <- drop1stloo(fit, scope$drop, training = training)
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
      nc <- match(c("loglik", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc], decreasing = TRUE)
      if (o[1L] == 1) 
        break
      change <- rownames(aod)[o[1L]]
    }
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)

    Terms <- terms(fit)
    bAIC <- extractAIC(fit)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    bloglik <- stloo.loglik(model = fit, training = training, nfit = fit)
    
    if (bloglik <= loglik - 1e-05) 
      break    

    nm <- nm + 1
    models[[nm]] <- list(fit = fit, deviance = mydeviance(fit), df.resid = n - 
                           edf, change = change, AIC = bAIC, loglik = bloglik)
  }
 
return(models[[nm]])  
}


drop1stloo <- function (object, scope, training) 
{
  tl <- attr(terms(object), "term.labels")
  if (missing(scope)) 
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L)) 
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>", 
                                                             scope), c("df", "loglik")))
  n <- length(object$residuals)
  edf <- n - object$df.residual
  ans[1, ] <- c(edf, stloo.loglik(object, training, object))
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq_len(ns)) {
    tt <- scope[i]
    
    nfit <- update(object, as.formula(paste("~ . -", tt)), 
                   evaluate = FALSE)
    ans[i + 1, ] <- c(edf + 1, stloo.loglik(object, training, nfit))
    
  }
  dfs <- ans[1L, 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, loglik = ans[, 2])
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
    dist.matrix <- as.matrix(x=dist(x=data[,c(lon.lab,lat.lab)], method = "euclidean",diag = T, upper = T))
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

