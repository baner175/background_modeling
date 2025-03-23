rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110
dat <- read.csv('Data/VBF_Cat2.csv', header = FALSE) # change data here

ni <- dat[,2]
N <- sum(ni)
k <- nrow(dat)
bins <- seq(l, u, length.out = k+1)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# Defining signal density and calculating signal region

fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

# Defining BW parameters
bkg_loc <- 91.2
bkg_scale <- 2.49/2

# generating pseudo unbinned data
set.seed(1234)
obs <- c()
for(i in 1:(length(bins))-1)
{
  obs <- c(obs, runif(ni[i], bins[i], bins[i+1]))
}

BW_power_likelihood <- function(power)
{
  total_mass <- integrate(f = function(x){
    1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-power))
  }, l, u)$value
  
  fi <- sapply(obs, function(x) {
    (1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-power)))/total_mass
  })
  return(-sum(log(fi)))
}
(qb_pow <- nlminb(start = 0.01,
                  objective = BW_power_likelihood,
                  lower = 0.005, upper = 10)$par)

qb_ <- function(x) 1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-qb_pow))
qb_mass <- integrate(qb_, l, u)$value
pow_mass <- integrate(function(x) x^(-qb_pow), l, u)$value
# BW_mass <- integrate(function(x) 1/((x-bkg_loc)^2 + bkg_scale^2), l, u)$value
qb <- function(x) (1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-qb_pow)))/qb_mass

obs_from_pow <- rtrunc(1e4, spec = 'pareto', a = l, b = u,
                       shape = qb_pow-1, scale = l)

obs_from_BW <- rtrunc(1e4, spec = 'cauchy', a = l, b = u,
                      location = bkg_loc, scale = bkg_scale)
bkg_only_data <- ifelse(rbinom(1e4, size = 1, prob = pow_mass/qb_mass),
                        obs_from_pow, obs_from_BW)

hist(obs, breaks = 50, probability = TRUE)
curve(dtrunc(x, spec = 'cauchy', a = l, b = u,
             location = bkg_loc, scale = bkg_scale), 
      l, u, add = TRUE, lwd = 2, col = 'red', lty = 1)
curve(dtrunc(x, spec = 'pareto', a = l, b = u,
             shape = qb_pow-1, scale = l), 
      l, u, add = TRUE, lwd = 2, col = 'blue', lty = 2)
curve(qb, l, u, add = TRUE, lwd = 2, col = 'purple', lty = 3)

# WITH THE BKG ONLY DATA FOR THE TIME BEING 
# WE ARE USING AN UNBINNED LIKELIHOOD

sf_bkg_likelihood <- function(eps, mu, data)
{
  fbi <- sapply(data, function(x){
    eps*dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mu, sd = sd_sig) + 
      (1-eps)*((1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-qb_pow)))/qb_mass)
  })
  cat(sprintf("\rEvaluating at mu: %f", mu))
  return(-sum(log(fbi)))
}

mu_seq <- seq(120, 130, length.out = 21)
eps_seq <- sapply(mu_seq, function(t)
{
  nlminb(start = 0.01,
         objective = sf_bkg_likelihood,
         lower = 0, upper = 1,
         mu = t, data = bkg_only_data)$par
})

(gb_sf_eps <- max(eps_seq))

sf_phys_likelihood <- function(eta, bin_ends, bin_counts)
{
  N_total <- sum(bin_counts)
  bi <- bin_ends
  mu_i <- c()
  for(i in 1:(length(bi)-1))
  {
    mu_i[i] <- N_total * integrate(f = function(x){
      eta*dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig) + 
        (1-eta)*(gb_sf_eps*dtrunc(x, spec = 'norm', a = l, b = u,
                                  mean = mean_sig, sd = sd_sig) + 
                   (1-gb_sf_eps)*((1/((x-bkg_loc)^2 + bkg_scale^2) +
                                     (x^(-qb_pow)))/qb_mass))
    }, lower = bi[i], upper = bi[i+1])$value
  }
  
  return(sum(mu_i - log(mu_i)*bin_counts))
}

safeguard_res <- nlminb(start = 0,
                        objective = sf_phys_likelihood,
                        lower = 0, upper = 1,
                        bin_ends = bins, bin_counts = ni)

# estimated proportion from safeguard:
(eta_sf <- safeguard_res$par)

ll0 <- -sf_phys_likelihood(eta = 0,bin_ends = bins,
                           bin_counts = ni)
ll1 <- -safeguard_res$objective

# safeguard p-value:
(p_val_sf <- 0.5*pchisq(-2*(ll0-ll1),
                        df = 1,
                        lower.tail = FALSE))

# estimated sample count and significance:
N_mid <- sum(ni[xi<=130 & xi>=120])
(S_hat <- eta_sf*N_mid)
(B_hat <- (1-eta_sf)*N_mid)

S_hat/sqrt(B_hat)
