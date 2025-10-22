library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)

#-------------------------------------------------------------------------------
## background fitting functions ##

qb_bkg_model_binned <- function(beta, nbins, bin_counts){
  k <- nbins
  bin_ends <- seq(l, u, length.out = k+1)
  qb_i <- sapply(1:k, function(i){
    integrate(function(x){
      dtrunc(x, spec = 'pareto', a = l, b = u,
             scale = l, shape = beta)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(bin_counts*log(qb_i)))
}

qb_bkg_model_unbinned <- function(beta, data){
  qb_i <- sapply(data, function(x){
    dtrunc(x, spec = 'pareto', a = l, b = u,
           scale = l, shape = beta)
  })
  return(-sum(log(qb_i)))
}

#-------------------------------------------------------------------------------
## simulation functions ##

simulated_power_binned <- function(eta, nbins, T_phys, lambda,
                                   mean1_in_gb, mean2_in_gb,
                                   sd_in_gb, nsims, seed = 12345,
                                   signif.level = 0.05,
                                   beta0 = NULL){
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f and lambda = %.2f', eta, lambda))
  cat("\n--------------------------------------\n")
  
  k <- nbins; B <- nsims
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[-1] + bin_ends[-(k+1)])/2
  
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig', 'bkg_rate', 'bkg_shape',
                      'qb_bkg_model_binned', 'qb_bkg_model_unbinned'))
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  if(is.null(beta0)){
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        set.seed(seeds[i])
        # total sample sizes:
        N <- rpois(1, lambda = T_phys)
        
        # physics-sample:
        s_samp <- rtrunc(N, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(N, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(N)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        ni <- sapply(1:k, function(i){
          sum((phys_samp>bin_ends[i])&(phys_samp<=bin_ends[i+1]))
        })
        
        opt <- nlminb(start = 0.01,
                      objective = qb_bkg_model_binned,
                      lower = 0, upper = Inf,
                      nbins = k, bin_counts = ni)
        beta_hat <- opt$par
        norm_S <- integrate(function(x) {
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          return(((fs/gb-1)^2)*gb)
        },l, u)$value |> sqrt()
        S2_vec <- sapply(xi,
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig)
                           qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta_hat)
                           fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean1_in_gb, sd = sd_in_gb)
                           fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean2_in_gb, sd = sd_in_gb)
                           gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
                           S_val <- fs/gb - 1
                           return((fs/gb-1)/(norm_S^2))
                         })
        theta_0_hat <- sum(S2_vec*ni)/N
        
        d_log_qb_xi <- sapply(xi, function(x){
          1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
        })
        d2_log_qb <- -1/(beta_hat^2) - 
          ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
          ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2
        
        d_normS2 <- -(1-2*lambda)*integrate(function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(((fs/gb)^2)*qb*d_log_qb)
        },l, u)$value
        
        d_S2 <- sapply(xi, function(x){ 
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          d_log_qb_xi <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          num <- -((norm_S^2)*(fs/(gb^2))*qb*(1-2*lambda)*d_log_qb_xi + (fs/gb-1)*d_normS2)
          denom <- norm_S^4
          return(num/denom)
        })
        
        J_hat <- -(1/k)*sum(ni*d2_log_qb)
        d_theta0_hat <- sum(d_S2*ni)/N
        var_S2_F_hat <- sum((S2_vec^2)*ni)/N - theta_0_hat^2
        c_hat <- N/k
        
        sig_theta0_hat <- sqrt(
          var_S2_F_hat + (d_theta0_hat^2)*((c_hat/J_hat)^2)*sum(ni*(d_log_qb_xi^2))/N + 
            2*d_theta0_hat*(c_hat/J_hat)*sum(ni*S2_vec*d_log_qb_xi)/N
        )
        
        sqrt(N)*theta_0_hat/sig_theta0_hat
      }
  }else{
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta0)
      fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                    mean = mean1_in_gb, sd = sd_in_gb)
      fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                    mean = mean2_in_gb, sd = sd_in_gb)
      gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
      S_val <- (fs/gb-1)
      return((S_val^2)*gb)
    },l, u)$value |> sqrt()
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        set.seed(seeds[i])
        # total sample sizes:
        N <- rpois(1, lambda = T_phys)
        
        # physics-sample:
        s_samp <- rtrunc(N, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(N, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(N)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        ni <- sapply(1:k, function(i){
          sum((phys_samp>bin_ends[i])&(phys_samp<=bin_ends[i+1]))
        })
        
        S2_vec <- sapply(xi,
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig)
                           qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta0)
                           fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean1_in_gb, sd = sd_in_gb)
                           fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                         mean = mean2_in_gb, sd = sd_in_gb)
                           gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
                           S_val <- fs/gb - 1
                           return(S_val/(norm_S^2))
                         })
        theta_0_hat <- sum(S2_vec*ni)/N
        
        sig_theta0_hat <- sqrt((1/N)*sum((S2_vec^2)*ni) - theta_0_hat^2)
        
        sqrt(N)*theta_0_hat/sig_theta0_hat
      }
  }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

simulated_power_unbinned <- function(eta, n_phys, lambda,
                                     mean1_in_gb, mean2_in_gb,
                                     sd_in_gb, nsims, seed = 12345,
                                     signif.level = 0.05,
                                     beta0 = NULL){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig', 'bkg_rate', 'bkg_shape',
                      'qb_bkg_model_binned', 'qb_bkg_model_unbinned'))
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  if(is.null(beta0)){
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        set.seed(seeds[i])
        
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
                      data = phys_samp)
        beta_hat <- opt$par
        norm_S <- integrate(function(x) {
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          return(((fs/gb-1)^2)*gb)
        },l, u)$value |> sqrt()
        d_normS2 <- -(1-2*lambda)*integrate(function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(((fs/gb)^2)*qb*d_log_qb)
        },l, u)$value
        
        S2_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                             mean = mean_sig, sd = sd_sig)
                                qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                             scale = l, shape = beta_hat)
                                fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                              mean = mean1_in_gb, sd = sd_in_gb)
                                fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                              mean = mean2_in_gb, sd = sd_in_gb)
                                gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
                                return((fs/gb-1)/(norm_S^2))
                              })
        d_S2 <- function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean1_in_gb, sd = sd_in_gb)
          fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                        mean = mean2_in_gb, sd = sd_in_gb)
          gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          num <- -((norm_S^2)*(fs/(gb^2))*qb*(1-2*lambda)*d_log_qb + (fs/gb-1)*d_normS2)
          denom <- norm_S^4
          return(num/denom)
        }
        
        d_log_qb_vec <- sapply(phys_samp, function(x){
          1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
        })
        
        theta_0_hat <- mean(S2_phys_vec)
        
        J_hat <- -(-1/(beta_hat^2) - 
                     ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
                     ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2)
        V_hat <- sum(d_log_qb_vec^2)/n_phys
        
        d_theta0_hat <- sapply(phys_samp, d_S2) |> mean()
        var_S2_F_hat <- mean(S2_phys_vec^2) - theta_0_hat^2
        sig_theta0_hat <- sqrt(
          var_S2_F_hat +  (d_theta0_hat^2)*V_hat/(J_hat^2) + 
            (2/J_hat)*d_theta0_hat*sum(S2_phys_vec*d_log_qb_vec)/n_phys
        )
        
        sqrt(n_phys)*theta_0_hat/sig_theta0_hat
      }
  }else{
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta0)
      fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                    mean = mean1_in_gb, sd = sd_in_gb)
      fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                    mean = mean2_in_gb, sd = sd_in_gb)
      gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
      S_val <- (fs/gb-1)
      return((S_val^2)*gb)
    },l, u)$value |> sqrt()
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        s_samp <- rtrunc(n_phys, a = l, b = u, spec = 'norm',
                         mean = mean_sig, sd = sd_sig)
        b_samp <- rtrunc(n_phys, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
        u_mask <- runif(n_phys)
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        S2_phys_vec <- sapply(phys_samp, 
                              function(x){
                                fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                             mean = mean_sig, sd = sd_sig)
                                qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                             scale = l, shape = beta0)
                                fs1 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                              mean = mean1_in_gb, sd = sd_in_gb)
                                fs2 <- dtrunc(x, a = l, b = u, spec = 'norm', 
                                              mean = mean2_in_gb, sd = sd_in_gb)
                                gb <- lambda*(fs1 + fs2) + (1-2*lambda)*qb
                                S_val <- fs/gb - 1
                                return((fs/gb-1)/(norm_S^2))
                              })
        theta_0_hat <- mean(S2_phys_vec)
        sig_theta0_hat<- sqrt(mean(S2_phys_vec^2) - theta_0_hat^2)
        
        sqrt(n_phys)*theta_0_hat/sig_theta0_hat
      }
  }
  stopCluster(cl)
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

