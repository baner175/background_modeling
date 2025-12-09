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

simulated_power_binned <- function(eta, nbins, T_phys,
                                   r, nsims, seed = 12345,
                                   signif.level = 0.05,
                                   beta0 = NULL){
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  k <- nbins; B <- nsims; T_bkg <- r*T_phys
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
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
        ni <- sapply(1:k, function(i){
          sum((phys_samp>bin_ends[i])&(phys_samp<=bin_ends[i+1]))
        })
        mi <- sapply(1:k, function(i){
          sum((bkg_samp>bin_ends[i])&(bkg_samp<=bin_ends[i+1]))
        })
        
        opt <- nlminb(start = 0.01,
                 objective = qb_bkg_model_binned,
                 lower = 0, upper = Inf,
                 nbins = k, bin_counts = mi)
        beta_hat <- opt$par
        norm_S <- integrate(function(x) {
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          return(((fs/qb-1)^2)*qb)
        },l, u)$value |> sqrt()
        S2_vec <- sapply(xi,
                         function(x){
                           fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig)
                           qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta_hat)
                           S_val <- fs/qb - 1
                           return((fs/qb-1)/(norm_S^2))
                         })
        theta_0_hat <- sum(S2_vec*ni)/N
        delta_0_hat <- sum(S2_vec*mi)/M
        
        d_log_qb_xi <- sapply(xi, function(x){
          1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
        })
        d2_log_qb <- -1/(beta_hat^2) - 
          ((log(l)^2)*l^(-beta_hat) - (log(u)^2)*u^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)) -
          ((log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat)))^2
        
        d_normS2 <- -integrate(function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          d_log_qb <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          return(((fs^2)/qb)*d_log_qb)
        },l, u)$value
        
        d_S2 <- sapply(xi, function(x){
          fs <- dtrunc(x, a = l, b = u, spec = 'norm', 
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta_hat)
          d_log_qb_xi <- 1/beta_hat - log(x) - (log(u)*u^(-beta_hat) - log(l)*l^(-beta_hat))/(l^(-beta_hat) - u^(-beta_hat))
          num <- -((norm_S^2)*(fs/qb)*d_log_qb_xi + (fs/qb-1)*d_normS2)
          denom <- norm_S^4
          return(num/denom)
        })
        
        J_hat <- -(1/k)*sum(mi*d2_log_qb)
        V_hat <- sum((d_log_qb_xi^2)*mi)/M
        cb_hat <- M/k
        d_theta0_hat <- sum(d_S2*ni)/N
        d_delta0_hat <- sum(d_S2*mi)/M
        d_theta0_T <- 1/(1-delta_0_hat)
        d_delta0_T <- (theta_0_hat-1)/((1-delta_0_hat)^2)
        cov_term <- (sum(mi*S2_vec*d_log_qb_xi)/M)
        
        var_S2_F_hat <- sum((S2_vec^2)*ni)/N - theta_0_hat^2
        var_S2_Fb_hat <- sum((S2_vec^2)*mi)/M - delta_0_hat^2
        
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
        
        test_num <- sqrt(M*N)*(eta_hat)
        
        denom1 <- (M/((1-delta_0_hat)^2)) * var_S2_F_hat
        
        denom2 <- N*(((theta_0_hat-1)^2)/((1-delta_0_hat)^4)) * var_S2_Fb_hat
        
        denom3 <- N*(V_hat/(J_hat^2))*(cb_hat^2)*((d_theta0_T*d_theta0_hat + d_delta0_T*d_delta0_hat)^2)
        
        denom4 <- (2*N*cb_hat/J_hat)*d_delta0_T*(d_theta0_T*d_theta0_hat + 
                                                   d_delta0_T*d_delta0_hat)*cov_term
        
        
        test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)
        
        test_num/test_denom
      }
  }else{
    norm_S <- integrate(function(x) {
      fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                   mean = mean_sig, sd = sd_sig)
      qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = beta0)
      S_val <- (fs/qb-1)
      return((S_val^2)*qb)
    },l, u)$value |> sqrt()
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
        phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
        
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
                           qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                                        scale = l, shape = beta0)
                           S_val <- fs/qb - 1
                           return((fs/qb-1)/(norm_S^2))
                         })
        theta_0_hat <- sum(S2_vec*ni)/N
        delta_0_hat <- sum(S2_vec*mi)/M
        
        sig_theta0_hat_sq <- (1/N)*sum((S2_vec^2)*ni) - theta_0_hat^2
        sig_delta0_hat_sq <- (1/M)*sum((S2_vec^2)*mi) - delta_0_hat^2
        
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
        
        test_num <- sqrt(N*M)*(eta_hat)
        test_denom <- sqrt(
          M*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
            N*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
        )
        test_num/test_denom
      }
  }
  
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

simulated_power_unbinned <- function(eta, n_phys,
                                   r, nsims, seed = 12345,
                                   signif.level = 0.05,
                                   beta0 = NULL){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
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
          eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
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
          
          
          test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
          
          denom1 <- n_bkg*var_S2_F_hat*(d_theta_T^2)
          denom2 <- n_phys*var_S2_Fb_hat*(d_delta_T^2)
          denom3 <- n_phys*(V_hat/J_hat^2)*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)^2
          denom4 <- 2*(n_phys/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta_hat + d_delta_T*d_delta_hat)
          
          test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)
          
          test_num/test_denom
      }
  }else{
    norm_S <- integrate(function(x) {
          fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                       mean = mean_sig, sd = sd_sig)
          qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
                       scale = l, shape = beta0)
          S_val <- (fs/qb-1)
          return((S_val^2)*qb)
        },l, u)$value |> sqrt()
    test_stat_eta <- foreach(i = 1:B, .combine = c,
                             .packages = c('truncdist', 'VGAM'),
                             .options.snow = opts) %dopar%
      {
        bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                           rate = bkg_rate, shape = bkg_shape)
        
        # physics-sample:
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
        
        eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
        
        test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
        test_denom <- sqrt(
          n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
            n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
        )
        test_num/test_denom
      }
  }
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

simulated_power_binned_uniform <- function(eta, nbins, T_phys,
                                   r, nsims, seed = 12345,
                                   signif.level = 0.05){
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  k <- nbins; B <- nsims; T_bkg <- r*T_phys
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[-1] + bin_ends[-(k+1)])/2
  
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig',
                      'qb_bkg_model_binned', 'qb_bkg_model_unbinned'))
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  norm_S <- integrate(function(x) {
    fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                 mean = mean_sig, sd = sd_sig)
    qb <- dunif(x, l, u)
    S_val <- (fs/qb-1)
    return((S_val^2)*qb)
  },l, u)$value |> sqrt()
  test_stat_eta <- foreach(i = 1:B, .combine = c,
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
      phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
      
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
                         qb <- dunif(x, l, u)
                         S_val <- fs/qb - 1
                         return((fs/qb-1)/(norm_S^2))
                       })
      theta_0_hat <- sum(S2_vec*ni)/N
      delta_0_hat <- sum(S2_vec*mi)/M
      
      sig_theta0_hat_sq <- (1/N)*sum((S2_vec^2)*ni) - theta_0_hat^2
      sig_delta0_hat_sq <- (1/M)*sum((S2_vec^2)*mi) - delta_0_hat^2
      
      eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
      
      test_num <- sqrt(N*M)*(eta_hat)
      test_denom <- sqrt(
        M*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
          N*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
      )
      test_num/test_denom
}
  
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

simulated_power_unbinned_uniform <- function(eta, n_phys,
                                     r, nsims, seed = 12345,
                                     signif.level = 0.05){
  
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  B <- nsims; n_bkg <- r*n_phys
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig',
                      'qb_bkg_model_binned', 'qb_bkg_model_unbinned'))
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  norm_S <- integrate(function(x) {
    fs <- dtrunc(x, a = l, b = u, spec = 'norm',
                 mean = mean_sig, sd = sd_sig)
    qb <- dunif(x, l, u)
    S_val <- (fs/qb-1)
    return((S_val^2)*qb)
  },l, u)$value |> sqrt()
  test_stat_eta <- foreach(i = 1:B, .combine = c,
                           .packages = c('truncdist', 'VGAM'),
                           .options.snow = opts) %dopar%
    {
      bkg_samp <- rtrunc(n_bkg, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
      
      # physics-sample:
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
      
      eta_hat <- (theta_0_hat - delta_0_hat)/(1-delta_0_hat)
      
      test_num <- sqrt(n_phys*n_bkg)*(eta_hat)
      test_denom <- sqrt(
        n_bkg*sig_theta0_hat_sq/((1- delta_0_hat)^2) + 
          n_phys*sig_delta0_hat_sq*((theta_0_hat-1)^2)/((1-delta_0_hat)^4)
      )
      test_num/test_denom
    }
  p_vals <- pnorm(test_stat_eta, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(simulated_power)
}

