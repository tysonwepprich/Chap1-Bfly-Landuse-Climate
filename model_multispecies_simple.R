# Is site suitable for a species?
# Observation error common for all species at a site x year
# eventually would like variables selection, SVSS or the like
# covariates (landuse and climate) need some work to minimize inputs into model
# can do model selection with WAIC to select fixed v random (not sure if this is automatic though)


source('dataprep4jags.R')

# Specify model in BUGS language
sink("MultiSpecies3.txt")
cat("
    model {
    for (i in 1:nsites){ # weird index because N and y different sizes due to Nt0
      for (k in 1:nspecies){
        # Latent suitability for site x species      
        z[i, k] ~ dbern(omega[k])

        for (j in 1:nyears){
          # Observation model
          # Zero-inflated poisson
          y[i, j, k] ~ dpois(N.eff[i, j, k] + 0.00001)
          N.eff[i, j, k] <- z[i, k] * N[i, j + 1, k]
    
          # Process model
          N[i, j + 1, k] <- N[i, j, k] * lam[i, j, k] # multiply growth rate by previous size
          log(lam[i, j, k]) <- beta0[k] + beta1[k] * log(N[i, j, k] + 1) + obsRE[i, j]

        } #nyears
      } #nspecies
    } #nsites
    
    # Priors
    # Suitability
    for (k in 1:nspecies){
    omega[k] ~ dunif(0, 1)
    beta0[k] ~ dnorm(mu.r.sp, tau.r.sp)
    beta1[k] ~ dnorm(mu.K.sp, tau.K.sp)
    }  

    # Observation RE
    for (i in 1:nsites){
      for (j in 1:nyears){
        obsRE[i, j] ~ dnorm(0, tau.obs)T(-10,10) #truncating helps convergence
      }
    }

    for (i in 1:nsites){
      for (k in 1:nspecies){
        N[i, 1, k] ~ dpois(Nstart[i, k])
      }
    }
    
    sigma.obs ~ dunif(0, 5)
    tau.obs <- pow(sigma.obs, -2)
    
    mu.r.sp ~ dnorm(0, 0.001)
    tau.r.sp <- pow(sigma.r.sp, -2)
    sigma.r.sp ~ dunif(0, 5)
    mu.K.sp ~ dnorm(0, 0.001)
    tau.K.sp <- pow(sigma.K.sp, -2)
    sigma.K.sp ~ dunif(0, 5)
    
    }
    
    ",fill = TRUE)
sink()


Ninit <- array(NA, dim = c(nsites, nyears + 1, nspecies))
zinit <- array(NA, dim = c(nsites, nspecies))
for (i in 1:nsites) {
  for (k in 1:nspecies) {
    zinit[i, k] <- ifelse(length(which(count_array[i, , k] > 0)) >= 1, 1, 0)
    Ninit[i, 1, k] <- rpois(1, mean(count_array[i, , k], na.rm = TRUE))
    for (j in 2:(nyears + 1)) {
      Ninit[i, j, k] <- NA
    }
  }
}



j.data <- list(y = count_array, nsites = nsites, nyears = nyears, nspecies = nspecies, Nstart = Nstartvals)

j.inits <- function(){
  list(
    z = zinit,
    # N = Ninit,
    omega = runif(nspecies, 0, 1),
    mu.K.sp = rnorm(1),
    mu.r.sp = rnorm(1),
    sigma.obs = runif(1, 0, 5), 
    sigma.r.sp = runif(1, 0, 5), 
    sigma.K.sp = runif(1, 0, 5)
  )
}

j.param <- c("omega", "mu.K.sp", "mu.r.sp", "sigma.obs", "sigma.r.sp", "sigma.K.sp", "N", "beta0", "beta1")

ni <- 100000
nt <- 150
nb <- 50000
nc <- 3

mod <- jags(j.data, inits = NULL, j.param, "MultiSpecies3.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
  
# library(devtools)
# install_github(repo='jagstools', username='johnbaums')
library(jagstools)

# show all results except for the many N nodes
round(jagsresults(x=mod, params='N', invert=TRUE), 2)
# 
# jags.Ns <- jagsresults(x=mod, param='N')
# N.ests <- array(data = round(jags.Ns[, "mean"]), dim = c(nsites, nyears, nspecies))
N.ests <- rearray(x = mod, param = "N", fields = "50%")
N.ests <- N.ests[[1]]
N.ests <- N.ests[, -1, ]

sp <- 1
plot(count_array[,,sp], N.ests[,,sp])
abline(0, 1)