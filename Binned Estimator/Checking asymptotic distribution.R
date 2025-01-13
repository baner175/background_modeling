rm(list = ls())

library(truncdist)
library(ggplot2)

############################################################################
###################### VERIFYING DELTA IS NEGATIVE #########################
############################################################################
l <- 1; u <- 5
mean_sig <- 2.5; sd_sig <- 0.1
mean_back <- 0.5; sd_back <- 2.5
rate_gb <- 0.35
eps <- 1e-3

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

find_d <- function(d)
{
  pl <- Fs(mean_sig-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(mean_sig+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

(M_lower <- mean_sig - r)
(M_upper <- mean_sig + r)

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 1.9*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;
gb <- function(x){
  fs_val_1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                     mean = mean1_in_gb, sd = sd_in_gb)
  fs_val_2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                     mean = mean2_in_gb, sd = sd_in_gb)
  gb_val <- dtrunc(x, spec = 'exp', a = l, b = u,
                   rate = rate_gb)
  
  return(0.02*fs_val_1 + 0.02*fs_val_2 + 0.96*gb_val)
}

integrate(gb, l, u)

fb_true <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                         mean = mean_back, sd = sd_back)

norm_S <- integrate(function(t) ((fs(t)/gb(t)-1)^2)*gb(t), l, u)$value |> sqrt()
S1 <- function(x)
{
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  return((fs/gb(x)-1)/norm_S)
}

(delta <- integrate(function(t) S1(t)*fb_true(t), l, u)$value)

curve(gb, from = l, to = u,n = 1000, col = 'red', lwd = 2, lty = 2,
      ylim = c(0.05,0.45))
curve(fb_true, l, u, col = 'brown', lwd = 2, lty = 3, add = TRUE)
curve(fs, l, u, col = 'skyblue', lwd = 2, lty = 1, add = TRUE)
abline(v = c(M_lower, M_upper), col = 'black', lwd = 2, lty = 2)



############################################################################
#################### SIMULATING THE NORMAL VARIABLE ########################
############################################################################

n_rep <- 1e4  # do it for 1e4
T_pois <- 5e2
k <- 100
n_samp <- rpois(1, lambda = T_pois)
eta_true <- 0
test_stat <- c()

t_now <- Sys.time()
for(i in 1:n_rep)
{
  sig_samp <- rtrunc(n_samp, spec = 'norm', a = l, b = u,
                     mean = mean_sig, sd = sd_sig)
  bkg_samp <- rtrunc(n_samp, spec = 'norm', a = l, b = u,
                     mean = mean_back, sd = sd_back)
  u_mask <- runif(n_samp)
  obs_new <- ifelse(u_mask<eta_true, sig_samp, bkg_samp)
  
  hs <- hist(obs_new, probability = TRUE, breaks = seq(l, u, length.out = k+1))
  ni <- hs$counts
  xi <- hs$mids
  # pi_hat <- ni/T_pois
  theta_hat_hat <- mean(S1(xi)*ni)
  test_stat[i] <- sqrt(k)*(theta_hat_hat - delta*T_pois/k)/sqrt(sum(ni*S1(xi)^2)/k)
  message(sprintf("Iteration: %d/%d", i, n_rep))
}

# Time elapsed:
Sys.time() - t_now

hist(test_stat, probability = TRUE, breaks = 50)
curve(dnorm, add = TRUE, lwd = 2)

qqnorm(test_stat)
qqline(test_stat)

