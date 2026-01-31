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
              help = 'Expected physics sample size', metavar = "number"),
  make_option(c('-r', '--bkg_phys'), type = "double", default = 1,
              help = 'bkg to physics sample size ratio', metavar = "number"),
  make_option(c('--mu_left'), type = "double", default = 1,
              help = 'left limit for SS search window', metavar = "number"),
  make_option(c('--mu_right'), type = "double", default = 2,
              help = 'left limit for SS search window', metavar = "number"),
  make_option(c('-b', '--beta'), type = "double", default = NULL,
              help = 'known value of the parameter in gb', metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

B <- as.numeric(opt$N_iter); n_samp <- as.numeric(opt$n_samp)
r <- as.numeric(opt$bkg_phys)
mu1 <- as.numeric(opt$mu_left); mu2 <- as.numeric(opt$mu_right)
eta_true <- as.numeric(opt$eta); beta0 <- as.numeric(opt$beta)


################################################################
################ SIGNAL AND SIGNAL REGION ######################
################################################################

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

q <- function(x, beta = beta0){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta)
}

mu_seq <- seq(mu1, mu2, length.out = 20)
neg_ll_eps <- function(eps, mu, data){
  fi <- sapply(data, function(t){
    eps*dtrunc(t, spec = 'norm', a = l, b = u,
               mean = mu, sd = sd_sig) + 
      (1-eps)*dtrunc(t, spec = 'pareto', a = l, b = u,
                     scale = l, shape = beta0)
  })
  return(-sum(log(fi)))
}

neg_ll <- function(eta, eps, data){
  fi <- sapply(data, function(t){
    eta*dtrunc(t, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig) + 
      (1-eta)*(
        eps*dtrunc(t, spec = 'norm', a = l, b = u,
                      mean = mean_sig, sd = sd_sig)
        + (1-eps)*dtrunc(t, spec = 'pareto', a = l, b = u,
                            scale = l, shape = beta0)
      )
  })
  return(-sum(log(fi)))
}


set.seed(1234)
seeds <- sample.int(.Machine$integer.max, B)
cl <- makeCluster(8)
registerDoSNOW(cl)
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

start_time <- Sys.time()
LRT_sim_res <- foreach(i = 1:B,
                      .packages = c('truncdist', 'VGAM'),
                      .options.snow = opts,
                      .combine = rbind) %dopar%
  {
    set.seed(seeds[i])
    
    bkg_samp <- rtrunc(r*n_samp, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
    s_samp <- rtrunc(n_samp, a = l, b = u, spec = 'norm',
                     mean = mean_sig, sd = sd_sig)
    b_samp <- rtrunc(n_samp, a = l, b = u, spec = 'gamma',
                     rate = bkg_rate, shape = bkg_shape)
    u_mask <- runif(n_samp)
    phys_samp <- ifelse(u_mask <= eta_true, s_samp, b_samp)
    
    eps_hat_seq <- sapply(mu_seq, function(m){
      nlminb(start = 0.01, 
             objective = neg_ll_eps,
             lower = 0, upper = 1,
             mu = m,
             data = bkg_samp)$par
    })
    eps_hat <- max(eps_hat_seq)
    
    gb_eps <- function(x){
      eps_hat*dtrunc(x, spec = 'norm', a = l, b = u,
                     mean = mean_sig, sd = sd_sig) + 
        (1-eps_hat)*dtrunc(x, spec = 'pareto', a = l, b = u,
                           scale = l, shape = beta0)
    }
    S <- function(x) {fs(x)/gb_eps(x) - 1}
    norm_S <- integrate(function(x) (S(x)^2)*gb_eps(x), l, u)$value |> sqrt()
    S1 <- function(x) S(x)/norm_S
    delta_val <- integrate(function(x) S1(x)*fb_true(x), l, u)$value
    eta_hat <- nlminb(start = 0.01,
                      objective = neg_ll,
                      lower = 0, upper = 1,
                      eps = eps_hat,
                      data = phys_samp)$par
    eta_hat_C <- eta_hat*(eta_hat>0)
    
    ll_eta_hat_C <- -neg_ll(eta = eta_hat_C, eps = eps_hat, 
                            data = phys_samp)
    ll_0 <- -neg_ll(eta = 0, eps = eps_hat, 
                    data = phys_samp)
    test_stat <- 2*(ll_eta_hat_C - ll_0)
    
    c(delta_val, eta_hat_C, test_stat)
  }
close(pb)

Sys.time() - start_time
stopCluster(cl)

LRT_sim_res <- as.data.frame(LRT_sim_res)
colnames(LRT_sim_res) <- c('delta', 'eta_hat', 'test_stat')

file_name <- paste0('/home/baner175/Desktop/background_modeling/simulations/',
                    'Results/LRT_SSmax',
                    '_B(', B,
                    ')_beta0(', beta0,
                    ')_n_samp(', n_samp,
                    ')_eta(', eta_true,
                    ')_mu_window_(', mu1, '-', mu2,
                    ').csv')
write.csv(LRT_sim_res, file_name,
          row.names = FALSE)