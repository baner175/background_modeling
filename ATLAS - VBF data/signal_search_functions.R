source('necessary_functions.R')

binned_signal_search_1 <- function(lambda, fs, qb,
                                   mu1, mu2, sd,
                                   search_region,
                                   count_region,
                                   bin_mids, bin_counts)
{
  xi <- bin_mids; ni <- bin_counts
  l <- search_region[1]; u <- search_region[2]
  
  N_mid <- sum(ni[xi<=count_region[2] & xi>=count_region[1]])
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
  eta_hat_binned <- theta_hat_binned/norm_S
  theta_check <- mean(S1_vec*ni)
  t_stat_theta <- k*theta_check/sqrt(sum(ni*S1_vec^2))
  p_val_binned <- pnorm(t_stat_theta, lower.tail = FALSE)
  S_hat <- N_mid*eta_hat_binned
  B_hat <- N_mid*(1-eta_hat_binned)
  signif <- S_hat/sqrt(B_hat)
  
  return(c(eta_hat_binned, p_val_binned, S_hat, signif))
}

binned_signal_search <- function(lambda_seq, ...)
{
  return(sapply(lambda_seq,
                function(t) binned_signal_search_1(t, ...)))
}

unbinned_signal_search_1 <- function(lambda, 
                                     fs, qb,
                                     mu1, mu2, sd,
                                     search_region,
                                     count_region, data)
{
  obs <- data
  n <- length(obs)
  l <- search_region[1]; u <- search_region[2]
  
  N_mid <- sum(obs<=count_region[2] & obs>=count_region[1])
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
  eta_hat <- theta/norm_S
  t_stat_theta <- theta/se_theta
  p_val <- pnorm(t_stat_theta, lower.tail = FALSE)
  S_hat <- N_mid*eta_hat
  B_hat <- N_mid*(1-eta_hat)
  signif <- S_hat/sqrt(B_hat)
  return(c(eta_hat, p_val, S_hat, signif))
}

unbinned_signal_search <- function(lambda_seq, ...)
{
  {
    return(sapply(lambda_seq,
                  function(t) unbinned_signal_search_1(t, ...)))
  }
}
