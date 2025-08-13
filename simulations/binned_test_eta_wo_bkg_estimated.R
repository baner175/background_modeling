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
              help = 'sample size', metavar = "number"),
  make_option(c('-l', '--lambda'), type = "double", default = 0,
              help = 'lambda in gb', metavar = "number"),
  make_option(c('--eps'),type = "double", default = 1e-3,
              help = '1 - signal region mass', metavar = "number"),
  make_option(c('-k', '--bins'), type = "integer", default = 100,
              help = 'Number of bins', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

eta_true <- as.numeric(opt$eta); B <- as.numeric(opt$N_iter)
T_phys <- as.numeric(opt$T_phys); lambda0 <- as.numeric(opt$lambda)
eps <- as.numeric(opt$eps); k <- as.numeric(opt$bins)

bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[-1] + bin_ends[-(k+1)])/2

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

qb <- function(x, beta){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb <- function(x, beta, lambda = lambda0) {
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

qb_bkg_model <- function(beta, bin_counts){
  qb_i <- sapply(1:k, function(i){
    integrate(function(x){
      dtrunc(x, spec = 'pareto', a = l, b = u,
             scale = l, shape = beta)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(bin_counts*log(qb_i)))
}

d_log_qb <- function(beta, x){
  1/beta - log(x) - (log(u)*u^(-beta) - log(l)*l^(-beta))/(l^(-beta) - u^(-beta))
}

beta_star <- uniroot(f = function(b) integrate(function(t) d_log_qb(b, t)*f(t), l, u)$value,
                     interval = c(1e-3,10))$root


norm_S_star <- integrate(function(x){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm',
               mean = mean_sig, sd = sd_sig)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta_star)
  gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
  S_val <- fs/gb - 1
  return((S_val^2)*gb)}, l, u)$value |> sqrt()


theta0_beta_star <- integrate(function(x){
  fs <- dtrunc(x, a = l, b = u, spec = 'norm',
               mean = mean_sig, sd = sd_sig)
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = beta_star)
  gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
  S_val <- fs/gb - 1
  S2_val <- S_val/(norm_S_star^2)
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
    
    N <- rpois(1, lambda = T_phys)
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
    
    opt <- tryCatch(
      nlminb(start = 0.01,
             objective = qb_bkg_model,
             lower = 0, upper = Inf,
             bin_counts = ni),
      error = function(e) return(NULL)
    )
    
    if (is.null(opt)) return(NA_real_)  # skip this iteration
    
    beta_hat <- opt$par
    
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
                   scale = l, shape = beta_hat)
      gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
      S_val <- fs/gb - 1
      return((S_val^2)*gb)
    }, l, u)$value |> sqrt()
    
    S2_vec <- sapply(xi, function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                        a = l, b = u,
                        spec = 'norm')
      fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                         a = l, b = u,
                         spec = 'norm')
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
      S_val <- (fs/gb - 1)
      return(S_val/(norm_S^2))
    })
    c_hat <- N/k
    d_log_qb_vec <- sapply(xi, function(x) d_log_qb(beta_hat, x))
    d2_log_qb <- -1/(beta_hat^2) - 
      ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
      ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2
    d_normS2 <- -(1-2*lambda0)*integrate(function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                        a = l, b = u,
                        spec = 'norm')
      fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                         a = l, b = u,
                         spec = 'norm')
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
      d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      return((((fs/gb)^2)*qb)*d_log_qb)
    },l, u)$value
    
    d_S2 <- sapply(xi, function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                        a = l, b = u,
                        spec = 'norm')
      fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                         a = l, b = u,
                         spec = 'norm')
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      gb <- lambda0*(fs_val1+fs_val2) + (1-2*lambda0)*qb
      d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      num <- -((norm_S^2)*(fs/(gb^2))*qb*d_log_qb*(1-2*lambda0) + (fs/gb-1)*d_normS2)
      denom <- norm_S^4
      return(num/denom)
    })
    
    J_hat <- -(1/k)*sum(ni*d2_log_qb)
    c_hat <- N/k
    d_theta_0_hat <- sum(d_S2*ni)/N
    theta0_hat <- sum(ni*S2_vec)/N

    test_num <- sqrt(N)*(theta0_hat - theta0_beta_star)
    
    test_denom <- sqrt(
      (sum((S2_vec^2)*ni)/N - theta0_hat^2) + 
        (d_theta_0_hat^2) * (c_hat^2) * (1/J_hat^2) * sum(ni*(d_log_qb_vec)^2)/N + 
        2*d_theta_0_hat*c_hat*(1/J_hat)*sum(ni*d_log_qb_vec*S2_vec)/N
    )
    
    test_num/test_denom
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/binned_test_eta_wo_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'lambda(',lambda0,')_',
                    'eta(',eta_true,')','.csv')

write.csv(data.frame(test_stat = test_stat),
          file_name, row.names = FALSE)