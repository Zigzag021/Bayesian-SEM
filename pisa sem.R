library(knitr)
library(lavaan)
library(blavaan)
library(tidyverse)
library(semPlot)
library(magrittr)
library(Matrix)
library(haven)
library(rethinking)
library(rjags)

# FULL SEM PISA DATA ----
pisa_data <- read.csv("enjoyment of reading.csv", header = TRUE)
complete_data <- pisa_data[which(
  pisa_data$ST160Q01IA != "NA" &
    pisa_data$ST160Q02IA != "NA" &
    pisa_data$ST160Q03IA != "NA" &
    pisa_data$ST160Q04IA != "NA" &
    pisa_data$ST160Q05IA != "NA" &
    pisa_data$PV1READ != "NA" &
    pisa_data$PV2READ != "NA" &
    pisa_data$PV3READ != "NA" &
    pisa_data$PV4READ != "NA" &
    pisa_data$PV5READ != "NA" &
    pisa_data$PV6READ != "NA" &
    pisa_data$PV7READ != "NA" &
    pisa_data$PV8READ != "NA" &
    pisa_data$PV9READ != "NA" &
    pisa_data$PV10READ != "NA" &
    pisa_data$ISCEDL != "NA"
)
, ]
complete_data$PV1 <- as.numeric(scale(complete_data$PV1READ))
complete_data$PV2 <- as.numeric(scale(complete_data$PV2READ))
complete_data$PV3 <- as.numeric(scale(complete_data$PV3READ))
complete_data$PV4 <- as.numeric(scale(complete_data$PV4READ))
complete_data$PV5 <- as.numeric(scale(complete_data$PV5READ))
complete_data$PV6 <- as.numeric(scale(complete_data$PV6READ))
complete_data$PV7 <- as.numeric(scale(complete_data$PV7READ))
complete_data$PV8 <- as.numeric(scale(complete_data$PV8READ))
complete_data$PV9 <- as.numeric(scale(complete_data$PV9READ))
complete_data$PV10 <- as.numeric(scale(complete_data$PV10READ))

sem_data <- list(
  N = nrow(complete_data),
  ST160Q01IA = complete_data$ST160Q01IA,
  ST160Q02IA = complete_data$ST160Q02IA,
  ST160Q03IA = complete_data$ST160Q03IA,
  ST160Q04IA = complete_data$ST160Q04IA,
  ST160Q05IA = complete_data$ST160Q05IA,
  PV1 = complete_data$PV1,
  PV2 = complete_data$PV2,
  PV3 = complete_data$PV3,
  PV4 = complete_data$PV4,
  PV5 = complete_data$PV5,
  ISCEDL = complete_data$ISCEDL
)


model_sem <- '
    enjoyment =~ ST160Q01IA + prior("normal(0.6,0.4)")*ST160Q02IA + prior("normal(0.6,0.4)")*ST160Q03IA + prior("normal(0.6,0.4)")*ST160Q04IA + prior("normal(0.6,0.4)")*ST160Q05IA
    reading =~ PV1 + prior("normal(0.6,0.4)")*PV2 + prior("normal(0.6,0.4)")*PV3 + prior("normal(0.6,0.4)")*PV4 + prior("normal(0.6,0.4)")*PV5
    reading ~ enjoyment + ISCEDL
'
semPaths(semPlotModel_lavaanModel(model_sem))

# jags model ----
sem_jags <- '
  model {
  for(i in 1:N) {
    ST160Q01IA[i] ~ dnorm(mu[i,1], 1/theta[1,1])
    ST160Q02IA[i] ~ dnorm(mu[i,2], 1/theta[2,2])
    ST160Q03IA[i] ~ dnorm(mu[i,3], 1/theta[3,3])
    ST160Q04IA[i] ~ dnorm(mu[i,4], 1/theta[4,4])
    ST160Q05IA[i] ~ dnorm(mu[i,5], 1/theta[5,5])
    PV1[i] ~ dnorm(mu[i,6], 1/theta[6,6])
    PV2[i] ~ dnorm(mu[i,7], 1/theta[7,7])
    PV3[i] ~ dnorm(mu[i,8], 1/theta[8,8])
    PV4[i] ~ dnorm(mu[i,9], 1/theta[9,9])
    PV5[i] ~ dnorm(mu[i,10], 1/theta[10,10])
    ST160Q02IAnew[i] ~ dnorm(mu[i,2], 1/theta[2,2])
    ST160Q03IAnew[i] ~ dnorm(mu[i,3], 1/theta[3,3])
    ST160Q04IAnew[i] ~ dnorm(mu[i,4], 1/theta[4,4])
    ST160Q05IAnew[i] ~ dnorm(mu[i,5], 1/theta[5,5])
    PV2new[i] ~ dnorm(mu[i,7], 1/theta[7,7])
    PV3new[i] ~ dnorm(mu[i,8], 1/theta[8,8])
    PV4new[i] ~ dnorm(mu[i,9], 1/theta[9,9])
    PV5new[i] ~ dnorm(mu[i,10], 1/theta[10,10])

    # lvs
    eta[i,1] ~ dnorm(mu_eta[i,1], 1/psi[1,1])
    eta[i,2] ~ dnorm(mu_eta[i,2], 1/psi[2,2])
  }

  # mu definitions
  for(i in 1:N) {
    mu[i,1] <- nu[1,1] + lambda[1,1]*eta[i,1]
    mu[i,2] <- nu[2,1] + lambda[2,1]*eta[i,1]
    mu[i,3] <- nu[3,1] + lambda[3,1]*eta[i,1]
    mu[i,4] <- nu[4,1] + lambda[4,1]*eta[i,1]
    mu[i,5] <- nu[5,1] + lambda[5,1]*eta[i,1]
    mu[i,6] <- nu[6,1] + lambda[6,2]*eta[i,2]
    mu[i,7] <- nu[7,1] + lambda[7,2]*eta[i,2]
    mu[i,8] <- nu[8,1] + lambda[8,2]*eta[i,2]
    mu[i,9] <- nu[9,1] + lambda[9,2]*eta[i,2]
    mu[i,10] <- nu[10,1] + lambda[10,2]*eta[i,2]

    mu_eta[i,1] <- alpha[1,1]
    mu_eta[i,2] <- alpha[2,1] + beta[2,1]*eta[i,1] + beta[2,3]*ISCEDL[i]
  }

  # Assignments from parameter vector & equality constraints
  lambda[1,1] <- 1
  lambda[2,1] <- parameter[1]
  lambda[3,1] <- parameter[2]
  lambda[4,1] <- parameter[3]
  lambda[5,1] <- parameter[4]
  lambda[6,2] <- 1
  lambda[7,2] <- parameter[5]
  lambda[8,2] <- parameter[6]
  lambda[9,2] <- parameter[7]
  lambda[10,2] <- parameter[8]
  beta[2,1] <- parameter[9]
  beta[2,3] <- parameter[10]
  theta[1,1] <- pow(parameter[11],-1)
  theta[2,2] <- pow(parameter[12],-1)
  theta[3,3] <- pow(parameter[13],-1)
  theta[4,4] <- pow(parameter[14],-1)
  theta[5,5] <- pow(parameter[15],-1)
  theta[6,6] <- pow(parameter[16],-1)
  theta[7,7] <- pow(parameter[17],-1)
  theta[8,8] <- pow(parameter[18],-1)
  theta[9,9] <- pow(parameter[19],-1)
  theta[10,10] <- pow(parameter[20],-1)
  psi[1,1] <- pow(parameter[21],-1)
  psi[2,2] <- pow(parameter[22],-1)
  psi[3,3] <- 0.0759233129007715
  nu[1,1] <- parameter[23]
  nu[2,1] <- parameter[24]
  nu[3,1] <- parameter[25]
  nu[4,1] <- parameter[26]
  nu[5,1] <- parameter[27]
  nu[6,1] <- parameter[28]
  nu[7,1] <- parameter[29]
  nu[8,1] <- parameter[30]
  nu[9,1] <- parameter[31]
  nu[10,1] <- parameter[32]
  alpha[3,1] <- 2.91722498379079
  alpha[1,1] <- 0
  alpha[2,1] <- 0

  # Priors
  parameter[1] ~ dnorm(0.6, 6.25)
  parameter[2] ~ dnorm(0.6, 6.25)
  parameter[3] ~ dnorm(0.6, 6.25)
  parameter[4] ~ dnorm(0.6, 6.25)
  parameter[5] ~ dnorm(0.6, 6.25)
  parameter[6] ~ dnorm(0.6, 6.25)
  parameter[7] ~ dnorm(0.6, 6.25)
  parameter[8] ~ dnorm(0.6, 6.25)
  parameter[9] ~ dnorm(0,1e-2)
  parameter[10] ~ dnorm(0,1e-2)
  parameter[11] ~ dgamma(1,.5)
  parameter[12] ~ dgamma(1,.5)
  parameter[13] ~ dgamma(1,.5)
  parameter[14] ~ dgamma(1,.5)
  parameter[15] ~ dgamma(1,.5)
  parameter[16] ~ dgamma(1,.5)
  parameter[17] ~ dgamma(1,.5)
  parameter[18] ~ dgamma(1,.5)
  parameter[19] ~ dgamma(1,.5)
  parameter[20] ~ dgamma(1,.5)
  parameter[21] ~ dgamma(1,.5)
  parameter[22] ~ dgamma(1,.5)
  parameter[23] ~ dnorm(0,1e-3)
  parameter[24] ~ dnorm(0,1e-3)
  parameter[25] ~ dnorm(0,1e-3)
  parameter[26] ~ dnorm(0,1e-3)
  parameter[27] ~ dnorm(0,1e-3)
  parameter[28] ~ dnorm(0,1e-3)
  parameter[29] ~ dnorm(0,1e-3)
  parameter[30] ~ dnorm(0,1e-3)
  parameter[31] ~ dnorm(0,1e-3)
  parameter[32] ~ dnorm(0,1e-3)}'

jags_sem <- jags.model(textConnection(sem_jags),
                       data = sem_data,
                       n.chains = 4)
sem.samps.pred <- jags.samples(
  jags_sem,
  variable.names = c(
    "ST160Q02IAnew",
    "ST160Q03IAnew",
    "ST160Q04IAnew",
    "ST160Q05IAnew",
    "PV2new",
    "PV3new",
    "PV4new",
    "PV5new"
  ),
  n.iter = 3000
)

update(jags_sem, n.iter = 1000)
samples.sem <- coda.samples(jags_sem,
                            variable.names = c("lambda", "beta"),
                            n.iter = 3000)
samps.sem <- as.data.frame(samples.sem[[1]])
precis(samps.sem, depth = 3, prob = .95)

# prior predictive check ----
# Number of simulations
nsim <- 500

# Define the prior distributions for the parameters
nu_1 <- rnorm(nsim, mean = 0, sd = sqrt(1000))
lambda <- rnorm(nsim, mean = 0.6, sd = sqrt(0.16))
beta1 <- rnorm(nsim, mean = 0, sd = sqrt(100))
beta2 <- rnorm(nsim, mean = 0, sd = sqrt(100))
psi_k <- rgamma(nsim, shape = 1, rate = 0.5)
theta <- rgamma(nsim, shape = 1, rate = 0.5)
z <- complete_data$ISCEDL
# Simulate eta from its prior
eta1 <- rnorm(nsim, mean = 0, sd = sqrt(psi_k))
eta2 <- sapply(1:nrow(complete_data), function(i) {
  rnorm(nsim,
        mean = (beta1 * eta1 + beta2 * z[i]),
        sd = sqrt(psi_k))
})

yprior <- sapply(1:nrow(complete_data), function(i) {
  rnorm(nsim,
        mean = (nu_1 + lambda * eta2),
        sd = sqrt(theta))
})

dens(
  yprior[1, ],
  ylim = c(0, 0.6),
  xlim = c(-100 , 100),
  adj = 1,
  col = "purple"
)
for (i in 2:500) {
  dens(
    yprior[i, ],
    ylim = c(0, 0.6),
    xlim = c(-100 , 100),
    adj = 1,
    col = "purple",
    add = T
  )
}

# MCMC ----
# Extract samples for mcmc check
samples.sem <- coda.samples(
  jags_sem,
  variable.names = c("eta", "mu", "lambda", "beta"),
  n.iter = 3000
)
selected_samples.sem <- samples.sem[, c("eta[1,1]", "lambda[2,1]" , "beta[2,1]", "mu[1,2]")]
# Plot trace plots for parameters
plot(selected_samples.sem, density = FALSE)

# Posterior predictive checks ----
# Set up the plotting area to display multiple plots
par(mfrow = c(2, 4))  # Adjust the layout to 2 rows and 4 columns
# List of variables to plot
variables <- c("PV2",
               "PV3",
               "PV4",
               "PV5",
               "ST160Q02IA",
               "ST160Q03IA",
               "ST160Q04IA",
               "ST160Q05IA")

# Loop through each variable to create the plots
for (var in variables) {
  # Plot the density for the actual data
  plot(
    density(complete_data[[var]], bw = 1),
    # Adjust the bandwidth
    adj = 1,
    col = "black",
    ylim = c(0, 0.5),
    main = paste(var)
  )
  
  # Overlay the predicted densities
  for (i in 1:1000) {
    lines(density(sem.samps.pred[[paste0(var, "new")]][, i, ]), col = "lightblue")
  }
  
  # Redraw the actual data density line
  lines(density(complete_data[[var]], bw = 1), col = "black") # Adjust the bandwidth
}

par(mfrow = c(1, 1)) #change layout back

# Fit indices ----
model_null <- '
    ST160Q01IA ~~ ST160Q01IA
    ST160Q02IA ~~ ST160Q02IA
    ST160Q03IA ~~ ST160Q03IA
    ST160Q04IA ~~ ST160Q04IA
    ST160Q05IA ~~ ST160Q05IA
    PV1 ~~ PV1
    PV2 ~~ PV2
    PV3 ~~ PV3
    PV4 ~~ PV4
    PV5 ~~ PV5
    ISCEDL ~~ ISCEDL
'
bfit_sem <- bsem(model_sem, data = complete_data, std.lv = F)
fit_null <- bsem(model_null, data = complete_data, std.lv = F)
fitindices <- blavFitIndices(bfit_sem, baseline.model = fit_null)

summary(fitindices, central.tendency = "mean", prob = .95)