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

# Defining q_b for safeguard
bkg_loc <- 91.2
bkg_scale <- 2.49/2

set.seed(123456)
bkg_only_data <- rtrunc(1e4, spec = 'cauchy', a = l, b = u,
                        location = bkg_loc, scale = bkg_scale)

# WITH THE BKG ONLY DATA FOR THE TIME BEING 
# WE ARE USING AN UNBINNED LIKELIHOOD
sf_bkg_likelihood <- function(eps, mu, data)
{
  fbi <- sapply(data, function(t){
    eps*dtrunc(t, spec = 'norm', a = l, b = u,
               mean = mu, sd = sd_sig) + 
      (1-eps)*dtrunc(t, spec = 'cauchy', a = l, b = u,
                     location = bkg_loc, scale = bkg_scale)
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
                   (1-gb_sf_eps)*dtrunc(x, spec = 'cauchy', a = l, b = u,
                                     location = bkg_loc, scale = bkg_scale))
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
