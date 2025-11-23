rm(list = ls())
library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)
library(optparse)
library(latex2exp)
l <- 1; u <- 2
#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02

eps <- 1e-3; beta <- 3.87

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
qb <- function(x){
  dtrunc(x, spec = 'pareto', a = l, b = u,
         scale = l, shape = beta0)
}

mean1_in_gb <- 0.5*M_lower + 0.5*mean_sig; sd_in_gb <- 4*sd_sig
mean2_in_gb <- 0.4*M_upper + 0.6*mean_sig

gb <- function(x, lambda) {
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

lambda_seq <- seq(0, 0.03, length.out = 4)
curve(fb_true, from = M_lower - 0.05, to = M_upper + 0.05, col = 'brown',
      lwd = 2, lty = 1, ylab = '')
abline(v = c(M_lower, M_upper), lty = 3, lwd = 2,
       col = 'black')

mycols <- c('black', 'red', 'green', 'purple')
for(i in 1:length(lambda_seq)){
  curve(gb(x, lambda = lambda_seq[i]),
        col = mycols[i],
        lwd = 2, lty = i+1, add = TRUE)
}
delta_seq <- sapply(lambda_seq, function(lam){
  normS <- integrate(function(x){
    fs <- fs(x); gb <- gb(x, lambda = lam)
    S <- (fs/gb - 1)
    return((S^2)*gb)
  }, l, u)$value |> sqrt()
  delta <- integrate(function(x){
    fs <- fs(x); gb <- gb(x, lambda = lam)
    S1 <- (fs/gb - 1)/normS
    return(S1*fb_true(x))
  }, l, u)$value
  return(delta)
})
legend('top',
       legend = c(
         TeX('$f_b(x)$'),
         TeX(
           paste0(
             sprintf('$g_b(\\lambda = %f);$', rev(lambda_seq)),
             sprintf('   $\\delta = %.4f$', rev(delta_seq))
           )
         )
       ),
       bty = 'n', lty = c(1, (length(lambda_seq)+1):2), 
       col = c('brown', rev(mycols)), lwd = 2)

