rm(list = ls())
## Global parameters ##

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

l <- 1; u <- 2

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)

gb <- function(x, beta = beta0){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}


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
              help = 'known value of beta in the background', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

B <- as.numeric(opt$N_iter); n_phys <- as.numeric(opt$n_phys)
bkg_to_phys_ratio <- as.numeric(opt$bkg_phys)
n_bkg <- n_phys*bkg_to_phys_ratio; beta0 <- as.numeric(opt$beta)


eta_seq <- c(0, 0.01, 0.03, 0.05)

# simulation results with known beta:

low_binned_powers <- sapply(eta_seq, function(e) {
  simulated_power_binned(eta = e, nbins = 1e2, T_phys = n_phys,
                         r = bkg_to_phys_ratio, nsims = B,
                         seed = 12345, signif.level = 0.05,
                         beta0 = beta0)
})

mid_binned_powers <- sapply(eta_seq, function(e) {
  simulated_power_binned(eta = e, nbins = 5e2, T_phys = n_phys,
                         r = bkg_to_phys_ratio, nsims = B,
                         seed = 12345, signif.level = 0.05,
                         beta0 = beta0)
})

high_binned_powers <- sapply(eta_seq, function(e) {
  simulated_power_binned(eta = e, nbins = 750, T_phys = n_phys,
                         r = bkg_to_phys_ratio, nsims = B,
                         seed = 12345, signif.level = 0.05,
                         beta0 = beta0)
})


unbinned_powers <- sapply(eta_seq, function(e) {
  simulated_power_unbinned(eta = e, n_phys = n_phys,
                           r = bkg_to_phys_ratio, nsims = B,
                           seed = 12345, signif.level = 0.05,
                           beta0 = beta0)
})

tab_beta_known <- data.frame(cbind(low_binned_powers, 
                                   mid_binned_powers,
                                   high_binned_powers,
                                   unbinned_powers))
colnames(tab_beta_known) <- paste0(c('low', 'mid', 'high', 'unbinned'), ' - beta known')

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/power_simulation_w_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,').csv')

write.csv(x = tab_beta_known, file_name,
          row.names = FALSE)