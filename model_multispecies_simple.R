# just species random effects for R intercept and K slopes. Site and year random effects on R intercept.

# Specify model in BUGS language
sink("MultiSpecies2.txt")
cat("
    model {
    for (i in 1:nsites){ # weird index because N and y different sizes due to Nt0
    for (k in 1:nspecies){
    # Latent suitability for site x species      
    z[i, k] ~ dbern(omega[k])
    
    for (j in 2:(nyears + 1)){
    # Observation model
    # Zero-inflated poisson
    y[i, j - 1, k] ~ dpois(N.eff[i, j - 1, k] + 0.00001)
    N.eff[i, j - 1, k] <- z[i, k] * N[i, j, k]
    
    # Process model
    N[i, j, k] <- N[i, j - 1, k] * lam[i, j - 1, k] # multiply growth rate by previous size
    log(lam[i, j - 1, k]) <- beta0[k] + beta1[k] * log(N[i, j - 1, k] + 1) + siteRE[i] + yearRE[j - 1] 
    #         + eps[i, j - 1, k]
    #         eps[i, j - 1, k] ~ dnorm(0, tau.proc[k])
    
    }
    }
    }
    
    # Priors
    # Suitability
    for (k in 1:nspecies){
    omega[k] ~ dunif(0, 1)
    beta0[k] ~ dnorm(mu.r.sp, tau.r.sp)
    beta1[k] ~ dnorm(mu.K.sp, tau.K.sp)
    }  
    
    for (i in 1:nsites){
    siteRE[i] ~ dnorm(0, tau.siteRE)
    for (k in 1:nspecies){
    N[i, 1, k] ~ dpois(Nstart[i, k])
    }
    }
    
    for (j in 1:nyears){
    yearRE[j] ~ dnorm(0, tau.yearRE)
    }
    
    #     for (k in 1:nspecies){
    #       sigma.proc[k] ~ dunif(0, 5)
    #       tau.proc[k] <- pow(sigma.proc[k], -2)
    #     }
    #     sigma.proc ~ dunif(0, 5)
    #     tau.proc <- pow(sigma.proc, -2)
    
    sigma.yearRE ~ dunif(0, 5)
    tau.yearRE <- pow(sigma.yearRE, -2)
    sigma.siteRE ~ dunif(0, 5)
    tau.siteRE <- pow(sigma.siteRE, -2)
    
    mu.r.sp ~ dnorm(0, 0.001)
    tau.r.sp <- pow(sigma.r.sp, -2)
    sigma.r.sp ~ dunif(0, 5)
    mu.K.sp ~ dnorm(0, 0.001)
    tau.K.sp <- pow(sigma.K.sp, -2)
    sigma.K.sp ~ dunif(0, 5)
    
    }
    
    ",fill = TRUE)
sink()


#add column of NAs for Nt0 in initial values of N
Ninit <- array(NA, dim = c(nsites, nyears))
# for (i in 1:nrow(Ninit)){
#   Ninit[i, ] <- rep(max(dat_array[i,,1], na.rm = TRUE) + 1, ncol(Ninit))
# }
Ninit <- cbind(rep(NA, nrow(Ninit)), Ninit)


j.data <- list(y = count_array, nsites = nsites, nyears = nyears, nspecies = nspecies, Nstart = Nstartvals)

j.inits <- function(){
  list(
    
    omega = runif(nspecies, 0, 1),
    mu.K.sp = rnorm(1),
    mu.r.sp = rnorm(1),
    sigma.yearRE = runif(1, 0, 5), 
    sigma.siteRE = runif(1, 0, 5), 
    sigma.r.sp = runif(1, 0, 5), 
    # sigma.proc = runif(nspecies, 0, 5),
    sigma.K.sp = runif(1, 0, 5)
  )
}
j.param <- c("omega", "mu.K.sp", "mu.r.sp", "sigma.yearRE", "sigma.siteRE", "sigma.r.sp", "N")

ni <- 10000
nt <- 8
nb <- 2000
nc <- 3

mod2MS <- jags(j.data, inits = NULL, j.param, "MultiSpecies2.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# library(devtools)
# install_github(repo='jagstools', username='johnbaums')
library(jagstools)

# show all results except for the many N nodes
round(jagsresults(x=mod2MS, params='N', invert=TRUE), 2)

jags.Ns <- jagsresults(x=mod2MS, param='N')
N.ests <- array(data = round(jags.Ns[, "mean"]), dim = c(nsites, nyears, nspecies))
rearray(x = mod2MS, param = "N", fields = "mean")

sp <- 8
plot(count_array[,,sp], N.ests[,,sp])
abline(0, 1)