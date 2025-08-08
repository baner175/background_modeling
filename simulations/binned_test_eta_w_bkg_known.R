rm(list = ls())

library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)
library(optparse)
l <- 1; u <- 2

option_list <- list(
  make_option(c("-e", "--eta"), type = "double", default = 0,
              help = "true value of eta", metavar = "number"),
  make_option(c("-B", "--N_iter"), type = "integer", default = 1e4,
              help = "Number of Iterations", metavar = "number"),
  make_option(c('-T', '--T_phys'), type = "integer", default = 2e3,
              help = 'Expected physics sample size', metavar = "number"),
  make_option(c('-r', '--bkg_phys'), type = "double", default = 1,
              help = 'bkg to physics sample size ratio', metavar = "number"),
  make_option(c('-k', '--bins'), type = "integer", default = 100,
              help = 'Number of bins', metavar = "number"),
  make_option(c('-b', '--beta'), type = "double", default = NULL,
              help = 'known value of beta in the background', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

eta_true <- as.numeric(opt$eta); B <- as.numeric(opt$N_iter);
k <- as.numeric(opt$bins)
T_phys <- as.numeric(opt$T_phys); bkg_to_phys_ratio <- as.numeric(opt$bkg_phys)
T_bkg <- T_phys*bkg_to_phys_ratio; beta0 <- as.numeric(opt$beta)

bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[-1] + bin_ends[-(k+1)])/2


#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
# true mixture density
f <- function(x) eta_true*fs(x)+(1-eta_true)*fb_true(x)

gb <- function(x, beta = beta0){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

norm_S <- integrate(function(x) {
  fs <- dtrunc(x, a = l, b = u, spec = 'norm',
               mean = mean_sig, sd = sd_sig)
  gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta0)
  S_val <- (fs/gb-1)
  return((S_val^2)*gb)
},l, u)$value |> sqrt()


set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)

cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
test_stat_eta <- foreach(i = 1:B, .combine = c,
                         .packages = c('truncdist', 'VGAM'),
                         .options.snow = opts) %dopar%
  {
    set.seed(seeds[i])
    # total sample sizes:
    M <- rpois(1, lambda = T_bkg); N <- rpois(1, lambda = T_phys)
    
    # bkg-only sample:
    bkg_samp <- rtrunc(M, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
    
    # physics-sample:
    s_samp <- rtrunc(N, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(N, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(N)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    ni <- sapply(1:k, function(i){
      sum((phys_samp>bin_ends[i])&(phys_samp<=bin_ends[i+1]))
    })
    mi <- sapply(1:k, function(i){
      sum((bkg_samp>bin_ends[i])&(bkg_samp<=bin_ends[i+1]))
    })
    
    S2_vec <- sapply(xi,
                     function(x){
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                         mean = mean_sig, sd = sd_sig)
                            gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                         scale = l, shape = beta0)
                            S_val <- fs/gb - 1
                            return((fs/gb-1)/(norm_S^2))
                          })
    theta_0_hat <- sum(S2_vec*ni)/N
    delta_0_hat <- sum(S2_vec*mi)/M
    
    sig_theta0_hat_sq <- (1/N)*sum((S2_vec^2)*ni) - theta_0_hat^2
    sig_delta0_hat_sq <- (1/M)*sum((S2_vec^2)*mi) - delta_0_hat^2
    
    eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    
    test_num <- sqrt(N*M)*(eta_hat - eta_true)
    test_denom <- sqrt(
      M*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
        N*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
    )
    test_num/test_denom
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/binned_test_eta_w_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_','eta(',eta_true,')','.csv')

write.csv(data.frame(test_stat = test_stat_eta),
          file_name, row.names = FALSE)
