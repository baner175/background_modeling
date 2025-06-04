source('necessary_functions.R')

binned_signal_search_1 <- function(lambda, fs, qb,
                                   mu1, mu2, sd,
                                   search_region,
                                   bin_mids, bin_counts)
{
  xi <- bin_mids; ni <- bin_counts
  l <- search_region[1]; u <- search_region[2]
  
  gb <- function(x){
    f1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu1, sd = sd)
    f2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu2, sd = sd)
    qb_val <- qb(x)
    return(lambda*(f1+f2) + (1-2*lambda)*qb_val)
  }
  
  norm_S <- integrate(function(t) {((fs(t)/gb(t) - 1)^2)*gb(t)},
                      l, u)$value |> sqrt()
  
  # Testing using the binned data:
  S1_vec <- sapply(xi, function(x) {
    f_sig <- fs(x)
    g_b <- gb(x)
    return((f_sig/g_b -1)/norm_S)
  })
  theta_hat_binned <- sum(S1_vec*ni)/sum(ni)
  eta_hat_binned <- (theta_hat_binned>=0)*(theta_hat_binned/norm_S)
  theta_check <- mean(S1_vec*ni)
  t_stat_theta <- k*theta_check/sqrt(sum(ni*S1_vec^2))
  p_val_binned <- pnorm(t_stat_theta, lower.tail = FALSE)
  S_hat <- sum(ni)*eta_hat_binned
  sigma_signif <- round(qnorm(p_val_binned, lower.tail = FALSE), 1)
  sigma_signif <- sigma_signif*(sigma_signif>0)
  
  return(c(eta_hat_binned,
           p_val_binned, 
           sigma_signif,
           round(S_hat, 2)))
}

binned_signal_search <- function(lambda_seq, ...)
{
  return(sapply(lambda_seq,
                function(t) binned_signal_search_1(t, ...)))
}


unbinned_signal_search_1 <- function(lambda, 
                                     fs, qb,
                                     mu1, mu2, sd,
                                     search_region, data)
{
  obs <- data
  n <- length(obs)
  l <- search_region[1]; u <- search_region[2]
  
  gb <- function(x){
    f1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu1, sd = sd)
    f2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu2, sd = sd)
    qb_val <- qb(x)
    return(lambda*(f1+f2) + (1-2*lambda)*qb_val)
  }
  norm_S <- integrate(function(t) {((fs(t)/gb(t) - 1)^2)*gb(t)},
                      l, u)$value |> sqrt()
  S1_vec <- sapply(obs, function(x) {
    f_sig <- fs(x)
    g_b <- gb(x)
    return((f_sig/g_b -1)/norm_S)
  })
  
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)
  eta_hat <- (theta>=0)*(theta/norm_S)
  t_stat_theta <- theta/se_theta
  p_val <- pnorm(t_stat_theta, lower.tail = FALSE)
  S_hat <- n*eta_hat
  sigma_signif <- round(qnorm(p_val, lower.tail = FALSE), 1)
  sigma_signif <- sigma_signif*(sigma_signif>0)
  return(c(eta_hat,
           p_val, 
           sigma_signif,
           round(S_hat, 2)))
}

unbinned_signal_search <- function(lambda_seq, ...)
{
  {
    return(sapply(lambda_seq,
                  function(t) unbinned_signal_search_1(t, ...)))
  }
}

unbinned_delta_estimation <- function(lambda, fs, qb,
                                      mu1, mu2, sd,
                                      search_region,
                                      bkg_data)
{
  obs <- bkg_data
  n <- length(obs)
  l <- search_region[1]; u <- search_region[2]
  
  gb <- function(x){
    f1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu1, sd = sd)
    f2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu2, sd = sd)
    qb_val <- qb(x)
    return(lambda*(f1+f2) + (1-2*lambda)*qb_val)
  }
  norm_S <- integrate(function(t) {((fs(t)/gb(t) - 1)^2)*gb(t)},
                      l, u)$value |> sqrt()
  S1_vec <- sapply(obs, function(x) {
    f_sig <- fs(x)
    g_b <- gb(x)
    return((f_sig/g_b -1)/norm_S)
  })
  
  delta_hat <- mean(S1_vec)
  se_delta <- sqrt((mean(S1_vec^2) - delta_hat^2)/n)
  cat(sprintf("\rEvaluating at lambda: %f", lambda))
  return(c(delta_hat, se_delta, pnorm(delta_hat/se_delta)))
}