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
  make_option(c('-n', '--n_samp'), type = "integer", default = 2e3,
              help = 'sample size', metavar = "number"),
  make_option(c('-l', '--lambda'), type = "double", default = 0,
              help = 'lambda in gb', metavar = "number"),
  make_option(c('-b', '--beta'), type = "double", default = NULL,
              help = 'known value of beta in the background', metavar = "number"),
  make_option(c('--eps'),type = "double", default = 1e-3,
              help = '1 - signal region mass', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

eta_true <- as.numeric(opt$eta); B <- as.numeric(opt$N_iter)
n_samp <- as.numeric(opt$n_samp); lambda0 <- as.numeric(opt$lambda)
beta0 <- as.numeric(opt$beta); eps <- as.numeric(opt$eps)

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)
Fs <- function(x) ptrunc(x, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig)

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))
r <- sol$root
M_lower <- mean_sig - r; M_upper <- mean_sig + r

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
# true mixture density
f <- function(x) eta_true*fs(x)+(1-eta_true)*fb_true(x)

qb <- function(x){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta0)
}

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb <- function(x) {
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta)
  
  return(lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb)
}


norm_S <- integrate(function(x){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm',
               mean = mean_sig, sd = sd_sig)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta0)
  gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
  S_val <- fs/gb - 1
  return((S_val^2)*gb)}, l, u)$value |> sqrt()


theta0 <- integrate(function(x){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm',
               mean = mean_sig, sd = sd_sig)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta0)
  gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
  S_val <- fs/gb - 1
  S2_val <- S_val/(norm_S^2)
  return(S2_val*f(x))
}, l, u)$value

set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)

cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
test_stat <- foreach(i = 1:B, .combine = c,
                     .packages = c('truncdist', 'VGAM'),
                     .options.snow = opts) %dopar%
  {
    set.seed(seeds[i])
    
    # physics-sample:
    s_samp <- rtrunc(n_samp, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_samp)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    
    S2_phys_vec <- sapply(phys_samp, 
                          function(x) {
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                         mean = mean_sig, sd = sd_sig)
                            fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                                              a = l, b = u,
                                              spec = 'norm')
                            fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                                               a = l, b = u,
                                               spec = 'norm')
                            qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                         scale = l, shape = beta0)
                            gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
                            S_val <- (fs/gb - 1)
                            return(S_val/(norm_S^2))
                          })
    
    theta0_hat <- mean(S2_phys_vec)
    
    test_num <- sqrt(n_samp)*(theta0_hat - theta0)
    test_denom <- sqrt(
      mean(S2_phys_vec^2) - theta0_hat^2
    )
    
    test_num/test_denom
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/unbinned_test_eta_wo_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'n_samp(',n_samp,')_',
                    'lambda(',lambda0,')_',
                    'eta(',eta_true,')','.csv')

write.csv(data.frame(test_stat = test_stat),
          file_name, row.names = FALSE)