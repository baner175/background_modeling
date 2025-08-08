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
  make_option(c('-n', '--n_phys'), type = "integer", default = 2e3,
              help = 'physics sample size', metavar = "number"),
  make_option(c('-r', '--bkg_phys'), type = "double", default = 1,
              help = 'bkg to physics sample size ratio', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

eta_true <- as.numeric(opt$eta); B <- as.numeric(opt$N_iter)
n_phys <- as.numeric(opt$n_phys); bkg_to_phys_ratio <- as.numeric(opt$bkg_phys)
n_bkg <- n_phys*bkg_to_phys_ratio


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

gb <- function(x, beta){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

gb_bkg_model <- function(beta, data){
  gb_i <- sapply(data, function(x){
    dtrunc(x, spec = 'pareto', a = l, b = u,
           scale = l, shape = beta)
  })
  return(-sum(log(gb_i)))
}

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
    # bkg-only sample:
    bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
    
    # physics-sample:
    s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_phys)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    opt <- tryCatch(
      nlminb(start = 0.01,
             objective = gb_bkg_model,
             lower = 0, upper = Inf,
             data = bkg_samp),
      error = function(e) return(NULL)
    )
    
    if (is.null(opt)) return(NA_real_)  # skip this iteration
    beta_hat <- opt$par
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      return(((fs/gb-1)^2)*gb)
    },l, u)$value |> sqrt()
    d_normS2 <- -integrate(function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      d_log_gb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      return(((fs^2)/gb)*d_log_gb)
    },l, u)$value
    
    S2_phys_vec <- sapply(phys_samp, 
                          function(x){
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean_sig, sd = sd_sig)
                            gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                         scale = l, shape = beta_hat)
                            return((fs/gb-1)/(norm_S^2))
                          })
    S2_bkg_vec <- sapply(bkg_samp,
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                        mean = mean_sig, sd = sd_sig)
                           gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta_hat)
                           return((fs/gb-1)/(norm_S^2))
                         })
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    J1_hat <- -(-1/(beta_hat^2) - 
      ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
      ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
    V1_hat <- sapply(bkg_samp,
                     function(x){
                       val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
                       return(val^2)
                     } ) |> mean()
    
    eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    
    d_theta_hat <- sapply(phys_samp, function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      S_val <- fs/gb - 1
      d_log_gb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      d_S <- (-fs/gb)*d_log_gb
      return(((norm_S^2)*d_S - S_val*d_normS2)/(norm_S^4))
    }) |> mean()
    d_delta_hat <- sapply(bkg_samp,function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      S_val <- fs/gb - 1
      d_log_gb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      d_S <- (-fs/gb)*d_log_gb
      return(((norm_S^2)*d_S - S_val*d_normS2)/(norm_S^4))
    }) |> mean()
    
    test_num <- sqrt(n_phys*n_bkg)*(eta_hat - eta_true)
    test_denom <- sqrt(
      (n_bkg/(1-delta_0_hat)^2) * (mean(S2_phys_vec^2) - theta_0_hat^2) +
        n_phys*((theta_0_hat - 1)^2/(1-delta_0_hat)^4) * (mean(S2_bkg_vec^2) - delta_0_hat^2) + 
        n_phys*V1_hat/(J1_hat^2) * (-d_theta_hat/(1-delta_0_hat) + ((theta_0_hat-1)/((1-delta_0_hat)^2))*d_delta_hat)^2 + 
        2 * (n_phys/J1_hat) * (theta_0_hat-1)/((1-delta_0_hat)^2) * 
        (-d_theta_hat/(1-delta_0_hat) + ((theta_0_hat-1)/((1-delta_0_hat)^2))*d_delta_hat) * 
        sapply(bkg_samp, function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          S2_val <- (fs/gb - 1)/(norm_S^2)
          d_log_gb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(S2_val*d_log_gb)
        }) |> mean())
    
    test_num/test_denom
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/unbinned_test_eta_w_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')

write.csv(data.frame(test_stat = test_stat_eta),
          file_name, row.names = FALSE)
