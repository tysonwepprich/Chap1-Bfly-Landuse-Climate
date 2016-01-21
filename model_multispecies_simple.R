# just species random effects for R intercept and K slopes. Site and year random effects on R intercept.

source('dataprep4jags.R')

# Specify model in BUGS language
sink("MultiSpecies3.txt")
cat("
    model {
    for (i in 1:nsites){ # weird index because N and y different sizes due to Nt0
      for (k in 1:nspecies){
      # Latent suitability for site x species      
      z[i, k] ~ dbern(omega[k])
    
        for (j in 2:(nyears + 1)){
        # Zero-inflated poisson with additional overdispersion
        y[i, j - 1, k] ~ dpois(N.eff[i, j - 1, k] * eDispersion[i, j - 1] + 0.00001)
        N.eff[i, j - 1, k] <- z[i, k] * N[i, j, k]

        # Process model
        N[i, j, k] <- N[i, j - 1, k] * lam[i, j - 1, k] # multiply growth rate by previous size
        log(lam[i, j - 1, k]) <- beta0[k] + beta1[i, k] * log(N[i, j - 1, k] + 1) + yearRE[j - 1] 
    
        }
      }
    }

    for (i in 1:nsites){
      for (j in 1:nyears){
        # Observation error shared between species at site x year for overdispersion
        eDispersion[i, j] ~ dgamma(1 / sDispersion^2, 1 / sDispersion^2)
      }
    }


    # Priors
    # Suitability
    for (k in 1:nspecies){
    omega[k] ~ dunif(0, 1)
    beta0[k] ~ dnorm(mu.r.sp, tau.r.sp)
    }  

    # Overdispersion
    sDispersion ~ dunif(0, 5)

    for (i in 1:nsites){
      for (k in 1:nspecies){
      N[i, 1, k] ~ dpois(Nstart[i, k])
      beta1[i, k] ~ dnorm(mu.K, tau.K)
      }
    }
    
    for (j in 1:nyears){
    yearRE[j] ~ dnorm(0, tau.yearRE)
    }
    sigma.yearRE ~ dunif(0, 5)
    tau.yearRE <- pow(sigma.yearRE, -2)

#     sigma.siteRE ~ dunif(0, 5)
#     tau.siteRE <- pow(sigma.siteRE, -2)
    
    mu.r.sp ~ dnorm(0, 0.001)
    tau.r.sp <- pow(sigma.r.sp, -2)
    sigma.r.sp ~ dunif(0, 5)

    mu.K ~ dnorm(0, 0.001)
    tau.K <- pow(sigma.K, -2)
    sigma.K ~ dunif(0, 5)
    
    }
    
    ",fill = TRUE)
sink()


#add column of NAs for Nt0 in initial values of N
# Ninit <- array(NA, dim = c(nsites, nyears))
# for (i in 1:nrow(Ninit)){
#   Ninit[i, ] <- rep(max(dat_array[i,,1], na.rm = TRUE) + 1, ncol(Ninit))
# }
# Ninit <- cbind(rep(NA, nrow(Ninit)), Ninit)


j.data <- list(y = count_array, nsites = nsites, nyears = nyears, nspecies = nspecies, Nstart = Nstartvals)

j.inits <- function(){
  list(
    sDispersion = runif(1, 0, 5),
    omega = runif(nspecies, 0, 1),
    mu.K = rnorm(1),
    mu.r.sp = rnorm(1),
    sigma.yearRE = runif(1, 0, 5), 
    sigma.r.sp = runif(1, 0, 5), 
    sigma.K = runif(1, 0, 5)
  )
}
j.param <- c("sDispersion", "omega", "mu.K", "mu.r.sp", "sigma.yearRE", 
             "sigma.r.sp", "sigma.K", "N")

ni <- 100000
nt <- 50
nb <- 50000
nc <- 3

mod2MS <- jags(j.data, inits = NULL, j.param, "MultiSpecies3.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# what about using rjags instead of R2jags?
library(rjags)
mod3MS <- jags.model("MultiSpecies3.txt", j.data, inits = NULL, n.chains = nc, n.adapt = 5000)

update( mod3MS , n.iter=nb )
codaSamples = coda.samples( mod3MS , variable.names = j.param,
                            n.iter=ni, thin = nt)
nChains <- 3
# parallel chains hopefully
library(runjags)
runJagsOut <- run.jags( method="rjags" ,
                        model="MultiSpecies3.txt" ,
                        monitor=j.param ,
                        data= j.data ,
                        inits= NULL ,
                        n.chains = 3,
                        adapt = 5000 ,
                        burnin = 10000 ,
                        sample= 20000 ,
                        thin=10 ,
                        summarise=TRUE ,
                        plots=TRUE )
codaSamples = as.mcmc.list( runJagsOut )


# show all results except for the many N nodes
round(jagsresults(x=mod2MS, params='N', invert=TRUE), 2)

jags.Ns <- jagsresults(x=mod2MS, param='N')
N.ests.2 <- array(data = round(jags.Ns[, "mean"]), dim = c(nsites, nyears + 1, nspecies))
test <- rearray(x = mod2MS, param = "N", fields = "mean")
N.ests <- round(test$N)

sp <- 3
plot(count_array[,,sp], N.ests[,-1,sp], xlim = c(0, 500), ylim = c(0, 500))
abline(0, 1)
