# Simulation to test different modeling approaches
set.seed(503)
nyr <- 20
nsite <- 50
nspec <- 30

# growth rate varies by species
r <- 0
rsig <- 0.01
rvec <- rnorm(nspec, r, rsig)

# from Chandler, det prob and iota for immigration into sites
iota <- 1
dp <- .7

# K varies by species and site
# Get starting values
Kmat <- array(data = NA, c(nsite, nspec))
Nst <- array(data = NA, c(nsite, nspec))
Kspec <- rpois(nspec, runif(nspec, 10, 500))

for (i in 1:nspec){
  for (j in 1:nsite){
    Ksite <- runif(nsite, min = .5, max = 2)
    Kmat[j, i] <- rpois(1, lambda = (Kspec[i] * Ksite[j])) + 1
    Nst[j, i] <- rpois(1, Kmat[j, i]) + 1 
  }
}

# weather effects vary by year and species
wspec <- rnorm(nspec, 0, .2)
wyr <- rnorm(nyr, mean = 0, sd = 1)

# simulate process with ricker model and weather effects and yearly noise
pops <- array(data = NA, dim = c(nsite, nyr, nspec))
counts <- array(data = NA, dim = c(nsite, nyr, nspec))
mu <- array(data = NA, dim = c(nsite, nyr, nspec))


# on log scale, having trouble with negative populations
# for (i in 1:nspec){
#   for (j in 1:nsite){
#     pops[j,1,i] <- log(Nst[j,i])
#     for (k in 2:nyr){
#       pops[j,k,i] <- pops[j,k-1,i] + (rvec[i] * (1 - (pops[j, k-1,i] / log(Kmat[j, i])))) + 
#                                            wspec[i] * wyr[k] + rnorm(1, 0, .1)
#     }
#   }
# }
# 
# for (i in 1:length(counts)){
#   counts[i] <- rpois(1, exp(pops[i]))
# }

# from Chandler, Ricker + immigration


for (i in 1:nspec){
  for (j in 1:nsite){
    pops[j,1,i] <- Nst[j,i]
    counts[j,1,i] <- rbinom(1, pops[j,1,i], dp)
    for (k in 2:nyr){
      mu[j,k,i] <- pops[j,k-1,i] * exp(rvec[i] * (1 - (pops[j, k-1,i] / Kmat[j, i])) + 
        wspec[i] * wyr[k] + rnorm(1, 0, .1)) + iota
      pops[j,k,i] <- rpois(1, mu[j,k,i])
      counts[j,k,i] <- rbinom(1, pops[j,k,i], dp)
    }
  }
}

dimnames(counts) <- list(paste("site", seq(1:nsite), sep = ""), 
                         seq(1:nyr), paste("spec", seq(1:nspec), sep = ""))
names(dimnames(counts)) <- list("site", "year", "species")

# different ways to analyze counts
library(reshape2)
library(data.table)
library(lme4)
test <- melt(counts, varnames = names(dimnames(counts)), value.name = "count")

cov <- data.frame(year = 1:nyr, cov = wyr)
dat <- merge(test, cov, by = "year")

# add last years count (slowly!)
Nt1 <- numeric(length = nrow(dat))
for (i in 1:nrow(dat)){
  temp <- dat$count[which(dat$site == dat$site[i] & dat$year == (dat$year[i] - 1) & dat$species == dat$species[i])]
  if (length(temp) == 0){
    Nt1[i] <- NA
  }else{
    Nt1[i] <- temp
  }
}
dat$Nt1 <- Nt1


# 1. model growth rates with complete cases
dat_edit <- data.table(dat)
dat_edit <- dat_edit[count > 0 | Nt1 > 0] #remove rows where both Nt and Nt-1 are zero.
dat_edit[, loglam := log(count + 1) - log(Nt1 + 1)]
dat_mod <- dat_edit[complete.cases(dat_edit)]
dat_mod[, logNt1 := log(Nt1 + 1)]

mod1 <- lmer(loglam ~ logNt1 + cov +
               (1 + logNt1 + cov|species), data = dat_mod)


mod2 <- lmer(loglam ~ logNt1 + cov + (0 + logNt1|site) +
               (1 + logNt1 + cov|species), data = dat_mod)

mod3 <- lmer(loglam ~  cov +
               (1 +  cov|species), data = dat_mod)

dat_mod[, Zsite_Nt1 := scale(logNt1), by = list(site, species)]
dat_mod[, Zcov := scale(cov)]

mod4 <- lmer(loglam ~ Zsite_Nt1 + Zcov +
               (1 + Zsite_Nt1 + Zcov|species), data = dat_mod)

library(MuMIn)
library(coefplot2)
coefplot2(mod1)

library(sjPlot)
sjp.setTheme(theme = theme_minimal())

sjt.lmer(mod1)
r.squaredGLMM(mod4)

sjp.lmer(mod1,facet.grid = FALSE,
         sort.coef = "sort.all")
sjp.lmer(mod1, "coef")

