rm(list = ls())
library(truncdist)
## Global parameters ##

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2

n_bins <- c(30, 50, 100)

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)

source('functions_for_power_simulation_w_bkg.R')

library(optparse)

option_list <- list(
  make_option(c("-B", "--N_iter"), type = "integer", default = 1e4,
              help = "Number of Iterations", metavar = "number"),
  make_option(c('-n', '--n_phys'), type = "integer", default = 2e3,
              help = 'physics sample size', metavar = "number"),
  make_option(c('-r', '--bkg_phys'), type = "double", default = 1,
              help = 'bkg to physics sample size ratio', metavar = "number"),
  make_option(c('-b', '--beta'), type = "double", default = NULL,
              help = 'known value of beta in the background', metavar = "number"),
  make_option(c('-e', '--eta'), type = "double", default = 0,
              help = 'proportion of signal', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

B <- as.numeric(opt$N_iter); n_phys <- as.numeric(opt$n_phys)
bkg_to_phys_ratio <- as.numeric(opt$bkg_phys); eta <- as.numeric(opt$eta)
n_bkg <- n_phys*bkg_to_phys_ratio; beta0 <- as.numeric(opt$beta)


# simulation results with known beta:

binned_res <- sapply(n_bins, function(k){
  simulated_power_binned(eta = eta, nbins = k, T_phys = n_phys,
                         r = bkg_to_phys_ratio, nsims = B,
                         seed = 12345, signif.level = 0.05,
                         beta0 = beta0)
})


unbinned_res <- simulated_power_unbinned(eta = eta, n_phys = n_phys,
                           r = bkg_to_phys_ratio, nsims = B,
                           seed = 12345, signif.level = 0.05,
                           beta0 = beta0)

power_res <- data.frame(t(c(binned_res, unbinned_res)))
colnames(power_res) <- paste0(c(paste0(n_bins, ' - bins'), 'unbinned'), ' (beta known)')

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/WBKG__',
                    'beta_known(',beta0,')_',
                    'eta(',eta,')_',
                    'n_bins(',paste0(n_bins, collapse = '-'),')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,').csv')

write.csv(x = power_res, file_name,
          row.names = FALSE)