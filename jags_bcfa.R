library(coda)
library(Matrix)
library(MASS)
library(rjags)

# Simulate data ----
J <- 1000
I <- 6
K <- 2
psi <- matrix(c(1, 0.5,
                0.5, 0.8), nrow = K)  
beta <- seq(1, 2, by = .2)

# loading matrix
Lambda <- cbind(c(1, 1.5, 2, 0, 0, 0), c(0, 0, 0, 1, 1.5, 2))

# error covariance
Theta <- diag(0.3, nrow = I)

# factor scores
eta <- mvrnorm(J, mu = c(0, 0), Sigma = psi)

# error term
epsilon <- mvrnorm(J, mu = rep(0, ncol(Theta)),Sigma = Theta)

dat <- tcrossprod(eta, Lambda) + epsilon
dat_cfa  <-  dat %>% as.data.frame() %>% setNames(c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6"))

jags_data <- list(
  N = nrow(dat_cfa),
  Y1 = dat_cfa$Y1,
  Y2 = dat_cfa$Y2,
  Y3 = dat_cfa$Y3,
  Y4 = dat_cfa$Y4,
  Y5 = dat_cfa$Y5,
  Y6 = dat_cfa$Y6,
  iden = diag(2)
)

# jags model ----
model_jags <- "
  model {
    for(i in 1:N) {
      Y1[i] ~ dnorm(mu[i,1], 1/theta[1,1])
      Y2[i] ~ dnorm(mu[i,2], 1/theta[2,2])
      Y3[i] ~ dnorm(mu[i,3], 1/theta[3,3])
      Y4[i] ~ dnorm(mu[i,4], 1/theta[4,4])
      Y5[i] ~ dnorm(mu[i,5], 1/theta[5,5])
      Y6[i] ~ dnorm(mu[i,6], 1/theta[6,6])
  
      # Latent variables bivariate normal distribution with mean 0 and precision matrix
      eta[i,1:2] ~ dmnorm(mu_eta[i,1:2], ibpsi[1:2,1:2])
    }
  
    # Mu definitions
    for(i in 1:N) {
      # mean of ov
      mu[i,1] <- nu[1,1] + lambda[1,1]*eta[i,1]
      mu[i,2] <- nu[2,1] + lambda[2,1]*eta[i,1]
      mu[i,3] <- nu[3,1] + lambda[3,1]*eta[i,1]
      mu[i,4] <- nu[4,1] + lambda[4,2]*eta[i,2]
      mu[i,5] <- nu[5,1] + lambda[5,2]*eta[i,2]
      mu[i,6] <- nu[6,1] + lambda[6,2]*eta[i,2]
      
      # mean of lv
      mu_eta[i,1] <- alpha[1,1]
      mu_eta[i,2] <- alpha[2,1]
    }
  
    # factor loading
    lambda[1,1] <- 1
    lambda[2,1] <- parameter[1]
    lambda[3,1] <- parameter[2]
    lambda[4,2] <- 1
    lambda[5,2] <- parameter[3]
    lambda[6,2] <- parameter[4]
    
    # variance of ov
    theta[1,1] <- pow(parameter[5], -1)
    theta[2,2] <- pow(parameter[6], -1)
    theta[3,3] <- pow(parameter[7], -1)
    theta[4,4] <- pow(parameter[8], -1)
    theta[5,5] <- pow(parameter[9], -1)
    theta[6,6] <- pow(parameter[10], -1)
    
    # intercept
    nu[1,1] <- parameter[11]
    nu[2,1] <- parameter[12]
    nu[3,1] <- parameter[13]
    nu[4,1] <- parameter[14]
    nu[5,1] <- parameter[15]
    nu[6,1] <- parameter[16]
    
    # mean of lv
    alpha[1,1] <- 0
    alpha[2,1] <- 0
    
    # variance covariance of lv
    psi[1,1] <- bpsi[1,1]
    psi[2,2] <- bpsi[2,2]
    psi[1,2] <- bpsi[1,2]
  
    # parameter priors
    parameter[1] ~ dnorm(0, 1e-2)
    parameter[2] ~ dnorm(0, 1e-2)
    parameter[3] ~ dnorm(0, 1e-2)
    parameter[4] ~ dnorm(0, 1e-2)
    parameter[5] ~ dgamma(1, .5)
    parameter[6] ~ dgamma(1, .5)
    parameter[7] ~ dgamma(1, .5)
    parameter[8] ~ dgamma(1, .5)
    parameter[9] ~ dgamma(1, .5)
    parameter[10] ~ dgamma(1, .5)
    parameter[11] ~ dnorm(0, 1e-3)
    parameter[12] ~ dnorm(0, 1e-3)
    parameter[13] ~ dnorm(0, 1e-3)
    parameter[14] ~ dnorm(0, 1e-3)
    parameter[15] ~ dnorm(0, 1e-3)
    parameter[16] ~ dnorm(0, 1e-3)
  
    for(k in 1:1) {
      # precision matrix Wishart distribution
      ibpsi[1:2,1:2] ~ dwish(iden, 3)
      
      # covariance matrix
      bpsi[1:2,1:2] <- inverse(ibpsi[1:2,1:2])
    }
  }"

# Compile jags model ----
jags_model <- jags.model(textConnection(model_jags), data = jags_data, n.chains = 3)
# Draw samples from the posterior
update(jags_model, n.iter = 1000)
samples <- coda.samples(jags_model, variable.names = c("lambda", "theta", "nu", "psi"), n.iter = 3000)
samps <- as.data.frame(samples[[1]])
precis(samps, depth = 3)
# Extract samples for lambda
lambda_samples <- samples[, c("lambda[1,1]", "lambda[2,1]", "lambda[3,1]", "lambda[4,2]", "lambda[5,2]", "lambda[6,2]")]
# Plot trace plots for lambda parameters
plot(lambda_samples, density = FALSE)

