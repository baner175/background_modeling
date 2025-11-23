rm(list = ls())
library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02
l <- 1; u <- 2
eta <- 0.02
B <- 1e3
n_phys <- 1e3; r <- 2
n_bkg <- r*n_phys
beta0 <- 3.87
#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5


qb_bkg_model_unbinned <- function(beta, data){
  qb_i <- sapply(data, function(x){
    dtrunc(x, spec = 'pareto', a = l, b = u,
           scale = l, shape = beta)
  })
  return(-sum(log(qb_i)))
}

set.seed(12345)
seeds <- sample.int(.Machine$integer.max, B)

cl <- makeCluster(8)
registerDoSNOW(cl)
clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig', 'bkg_rate', 'bkg_shape',
                    'qb_bkg_model_unbinned'))
pb <- txtProgressBar(max = B, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

eta_res <- foreach(i = 1:B, .combine = rbind,
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
    phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
    
    opt <- nlminb(start = 0.01,
                  objective = qb_bkg_model_unbinned,
                  lower = 0, upper = Inf,
                  data = bkg_samp)
    beta_hat <- opt$par
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      return(((fs/qb-1)^2)*qb)
    },l, u)$value |> sqrt()
    d_normS2 <- -integrate(function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      return(((fs^2)/qb)*d_log_qb)
    },l, u)$value
    
    S2_phys_vec <- sapply(phys_samp, 
                          function(x){
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean_sig, sd = sd_sig)
                            qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                         scale = l, shape = beta_hat)
                            return((fs/qb-1)/(norm_S^2))
                          })
    S2_bkg_vec <- sapply(bkg_samp,
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                        mean = mean_sig, sd = sd_sig)
                           qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta_hat)
                           return((fs/qb-1)/(norm_S^2))
                         })
    d_S2 <- function(x){
      fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta_hat)
      d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      num <- -((norm_S^2)*(fs/qb)*d_log_qb + (fs/qb-1)*d_normS2)
      denom <- norm_S^4
      return(num/denom)
    }
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    J_hat <- -(-1/(beta_hat^2) - 
                 ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
                 ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
    V_hat <- sapply(bkg_samp,
                    function(x){
                      val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
                      return(val^2)
                    } ) |> mean()
    eta_hat_est <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    d_theta_hat <- sapply(phys_samp, d_S2) |> mean()
    d_delta_hat <- sapply(bkg_samp, d_S2) |> mean()
    d_theta_T <- 1/(1-delta_0_hat)
    d_delta_T <- (theta_0_hat - 1)/((1-delta_0_hat)^2)
    d_log_qb_bkg <- sapply(bkg_samp, function(x){
      val <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
      return(val)
    })
    cov_term <- mean(S2_bkg_vec*d_log_qb_bkg)
    var_S2_F_hat <- mean(S2_phys_vec^2) - theta_0_hat^2
    var_S2_Fb_hat <- mean(S2_bkg_vec^2) - delta_0_hat^2
    denom1 <- n_bkg*var_S2_F_hat*(d_theta_T^2)
    denom2 <- n_phys*var_S2_Fb_hat*(d_delta_T^2)
    denom3 <- n_phys*(V_hat/J_hat^2)*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)^2
    denom4 <- 2*(n_phys/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)
    test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)
    std_err_est <- test_denom/sqrt(n_phys + n_bkg)
    
    
    
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta0)
      S_val <- (fs/qb-1)
      return((S_val^2)*qb)
    },l, u)$value |> sqrt()
    S2_phys_vec <- sapply(phys_samp, 
                          function(x){
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                         mean = mean_sig, sd = sd_sig)
                            gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                         scale = l, shape = beta0)
                            S_val <- fs/gb - 1
                            return((fs/gb-1)/(norm_S^2))
                          })
    S2_bkg_vec <- sapply(bkg_samp, 
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig)
                           gb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta0)
                           S_val <- fs/gb - 1
                           return((fs/gb-1)/(norm_S^2))
                         })
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    sig_theta0_hat_sq <- mean(S2_phys_vec^2) - theta_0_hat^2
    sig_delta0_hat_sq <- mean(S2_bkg_vec^2) - delta_0_hat^2
    eta_hat_knw <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    test_denom <- sqrt(
      n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
        n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
    )
    std_err_knw <- test_denom/sqrt(n_phys+n_bkg)
    
    
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dunif(x, l, u)
      S_val <- (fs/qb-1)
      return((S_val^2)*qb)
    },l, u)$value |> sqrt()
    S2_phys_vec <- sapply(phys_samp, 
                          function(x){
                            fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                         mean = mean_sig, sd = sd_sig)
                            gb <- dunif(x, l, u)
                            S_val <- fs/gb - 1
                            return((fs/gb-1)/(norm_S^2))
                          })
    S2_bkg_vec <- sapply(bkg_samp, 
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig)
                           gb <- dunif(x, l, u)
                           S_val <- fs/gb - 1
                           return((fs/gb-1)/(norm_S^2))
                         })
    theta_0_hat <- mean(S2_phys_vec)
    delta_0_hat <- mean(S2_bkg_vec)
    sig_theta0_hat_sq <- mean(S2_phys_vec^2) - theta_0_hat^2
    sig_delta0_hat_sq <- mean(S2_bkg_vec^2) - delta_0_hat^2
    eta_hat_missp <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
    test_denom <- sqrt(
      n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
        n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
    )
    std_err_missp <- test_denom/sqrt(n_phys+n_bkg)
    
    c(eta_hat_est, eta_hat_knw, eta_hat_missp, std_err_est, std_err_knw, std_err_missp)
  }

eta_res <- as.data.frame(eta_res)
colnames(eta_res) <- c('eta_hat_est',
                       'eta_hat_knw',
                       'eta_hat_missp',
                       'se_est',
                       'se_knw',
                       'se_missp')

write.csv(eta_res, 'Results/eta_hat_se_comparison.csv',
          row.names = FALSE)

eta_res <- read.csv('Results/eta_hat_se_comparison.csv',
                    header = TRUE)
View(eta_res)
