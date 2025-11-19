rm(list = ls())
library(truncdist)
## Global parameters ##

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2

n_bins <- c(30, 50, 100)

# signal density and CDF
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig)

##parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)

source('functions_for_power_simulation_wobkg.R')

library(optparse)

option_list <- list(
  make_option(c("-B", "--N_iter"), type = "integer", default = 1e4,
              help = "Number of Iterations", metavar = "number"),
  make_option(c('-n', '--n_phys'), type = "integer", default = 2e3,
              help = 'physics sample size', metavar = "number"),
  make_option(c('-l', '--lambda'), type = "double", default = 0,
              help = 'lambda in gb', metavar = "number"),
  make_option(c('-b', '--beta'), type = "double", default = NULL,
              help = 'known value of beta in the background', metavar = "number"),
  make_option(c('-e', '--eta'), type = "double", default = 0,
              help = 'proportion of signal', metavar = "number"),
  make_option('--eps', type = "double", default = 1e-3,
              help = 'mass outside of the signal region', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

B <- as.numeric(opt$N_iter); n_phys <- as.numeric(opt$n_phys)
eta <- as.numeric(opt$eta); lambda <- as.numeric(opt$lambda)
beta0 <- as.numeric(opt$beta); eps <- as.numeric(opt$eps)

# Calculating signal region:

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig - d)
  pu <- Fs(mean_sig + d)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

M_lower <- mean_sig - r
M_upper <- mean_sig + r

# PROPOSAL BACKGROUND DENSITY PARAMETERS:
mean1_in_gb <- 0.5*M_lower + 0.5*mean_sig
mean2_in_gb <- 0.4*M_upper + 0.6*mean_sig
sd_in_gb <- 4*sd_sig

# simulation results with known beta:

binned_res <- sapply(n_bins, function(k){
  simulated_power_binned(eta = eta, nbins = k, 
                         T_phys = n_phys, 
                         lambda = lambda,
                         mean1_in_gb = mean1_in_gb, 
                         mean2_in_gb = mean2_in_gb,
                         sd_in_gb = sd_in_gb, 
                         nsims = B, seed = 12345,
                         signif.level = 0.05,
                         beta0 = beta0)
})

unbinned_res <- simulated_power_unbinned(eta = eta, n_phys = n_phys, 
                                         lambda = lambda,
                                         mean1_in_gb = mean1_in_gb, 
                                         mean2_in_gb = mean2_in_gb,
                                         sd_in_gb = sd_in_gb, 
                                         nsims = B, seed = 12345,
                                         signif.level = 0.05,
                                         beta0 = beta0)

normS <- integrate(function(x){
  fs <- fs(x)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta0)
  gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
  S <- (fs/gb - 1)
  return((S^2)*gb)
}, l, u)$value |> sqrt()
delta_val <- integrate(function(x){
  fs <- fs(x)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta0)
  gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
  S1 <- (fs/gb-1)/normS
  return(S1*fb_true(x))
}, l, u)$value

power_res <- data.frame(t(c(binned_res, unbinned_res, delta_val)))
colnames(power_res) <- c(paste0(c(paste0(n_bins, ' - bins'), 'unbinned'), ' (beta known)'),
                         'delta')

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/WOBKG__',
                    'beta_known(',beta0,')_',
                    'eta(',eta,')_',
                    'lambda(', lambda,')_',
                    'n_bins(',paste0(n_bins, collapse = '-'),')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,').csv')

write.csv(x = power_res, file_name,
          row.names = FALSE)