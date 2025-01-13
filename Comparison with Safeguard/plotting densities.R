rm(list = ls())

library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

################################################################
######################### SEARCH REGION ########################
################################################################

l <- 1; u <- 2

################################################################
################ SIGNAL AND SIGNAL REGION ######################
################################################################

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02
eps <- 1e-3

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)
fs <- Vectorize(fs)
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

(M_lower <- mean_sig - r)
(M_upper <- mean_sig + r)

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps

################################################################
####################### TRUE BACKGROUND ########################
################################################################

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)
#sample size
n <- 5e3

set.seed(12344)
#generating data from true background
bkg_data <- rtrunc(n, a = l, b = u, spec = 'gamma', rate = bkg_rate, shape = bkg_shape)


################################################################
############### CONSTRUCTING SAFEGUARD BACKGROUND ##############
################################################################

qb_sf_mod <- function(pars, data)
{
  beta <- pars
  qb_sf <- dtrunc(data, spec = 'pareto', a = l, b = u,
                  shape = beta, scale = l)
  return(-sum(log(qb_sf)))
}

(qb_beta <- nlminb(start = 1,
                   objective = qb_sf_mod,
                   lower = 0, upper = Inf,
                   data = bkg_data)$par) # 3.867429

hist(bkg_data, probability = TRUE, breaks = 70)
curve(fb_true, lwd = 2.5, col = 'red', add = TRUE, lty = 2)
curve(dtrunc(x, spec = 'pareto', a = l, b = u,
             shape = qb_beta, scale = l),
      lwd = 2.5, col = 'brown', add = TRUE)

qb_sf_sim_data <- rtrunc(n, a = l, b = u, spec = 'pareto',
                         scale = l, shape = qb_beta)

# guess for the smooth bkg
qb_sf <- function(x) dtrunc(x,a=l,b=u,spec="pareto",scale=l,shape = qb_beta)
qb_sf <- Vectorize(qb_sf)

bkg_mod <- function(pars, data, mean = mean_sig)
{
  lam <- pars
  fi <- lam*dtrunc(data, a = l, b = u, spec = 'norm',
                   mean = mean, sd = sd_sig) + 
    (1-lam)*dtrunc(data,a=l,b=u,spec="pareto",scale=l,shape=qb_beta)
  return(-sum(log(fi)))
}

mean_seq <- c(seq(M_lower, mean_sig, length.out = 25),
              seq(mean_sig, M_upper, length.out = 25)[-1])

res_sf_list <- sapply(mean_seq, function(x) nlminb(start = 0.1, objective = bkg_mod,
                                                   lower = 0, upper = 1,
                                                   data = qb_sf_sim_data, mean = x)$par)

(lam_sf_hat <- max(res_sf_list)) # 0.002629639

gb_sf <- function(x) {lam_sf_hat*dtrunc(x, a = l, b = u, spec = 'norm',
                                        mean = mean_sig, sd = sd_sig) + 
    (1-lam_sf_hat)*dtrunc(x,a=l,b=u,spec="pareto",scale=l,shape=qb_beta)}
gb_sf <- Vectorize(gb_sf)

# This set up 
S_sf <- function(x) {fs(x)/gb_sf(x) - 1}
norm_S_sf <- integrate(function(x) S_sf(x)*S_sf(x)*gb_sf(x), l, u)$value |> sqrt()
S1_sf <- function(x) S_sf(x)/norm_S_sf
(delta_sf <- integrate(function(x) S1_sf(x)*fb_true(x), l, u)$value)


hist(bkg_data, probability = TRUE, breaks = 70)
curve(gb_sf, col = 'cyan', lwd = 2.5, from = l, to = u, add = TRUE)
curve(fs, col = 'blue', lwd = 2.5, from = l, to = u, add = TRUE)
curve(fb_true, col = 'brown', lwd = 2.5, from = l, to = u, add = TRUE)
legend('topright', col = c('cyan', 'brown', 'blue'),
       legend = c('proposed background from safeguard', 'true background', 'signal'),
       bty = 'n', lty = 1, lwd = 2.5)
abline(v = c(M_lower, M_upper), lwd = 2.5, col = 'black')

# negative log-likelihood for safeguard:
mix_mod_sf <- function(pars, data)
{
  eta <- pars
  mix_dens <- log(eta*dtrunc(data, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig) +
                    (1-eta)*(lam_sf_hat*dtrunc(data, a = l, b = u, spec = 'norm',
                                               mean = mean_sig, sd = sd_sig) + 
                               (1-lam_sf_hat)*dtrunc(data,a=l,b=u,spec="pareto",
                                                     scale=l,shape=qb_beta)))
  log_lik <- sum(mix_dens)
  return(-log_lik)
}

################################################################
###################### OUR METHOD ##############################
################################################################

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    a = l, b = u,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     a = l, b = u,
                     spec = 'norm')
  qb <- dtrunc(x, spec = 'pareto', a = l, b = u,
               scale = l, shape = qb_beta)
  
  return(fs_prop*(fs_val1+fs_val2) + (1-2*fs_prop)*qb)
}

plot(x = seq(1.2,1.4,length.out = 50),
     y = seq(0.5,2.5,length.out = 50), type = 'n',
     ylab = 'Proposal background densities',
     xlab = 'x')
curve(gb_test(x, fs_prop = 0), add = TRUE,
      col = 'black', lty = 1,
      lwd = 2.5,
      from = l, to = u, ylab = 'gb_test')

curve(gb_test(x, fs_prop = 0.002),
      col = 'blue', lty = 2,
      lwd = 2.5,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.005),
      col = 'green', lty = 3,
      lwd = 2.5,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.007),
      col = 'red', lty = 4,
      lwd = 2.5,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.01),
      col = 'orange', lty = 5,
      lwd = 2.5,
      add = TRUE)
curve(gb_sf,
      col = 'cyan', lty = 6,
      lwd = 2.5, add = TRUE)
curve(fs,
      col = 'skyblue', lty = 9,
      lwd = 2.5, add = TRUE)
curve(fb_true,
      col = 'brown', lty = 10,
      lwd = 2.5, add = TRUE)
abline(v = c(M_upper, M_lower), lwd = 2.5, col = 'black')


legend('topright', col = c('cyan', 'black', 'blue', 'green',
                           'red', 'orange', 'brown', 'skyblue'),
       lty = c(6,1:5,10,9), bty = 'n', lwd =2.5,
       legend=c('safeguard',
                TeX(sprintf(r'($g_b(\lambda = %f)$)',
                            c(0, 0.002, 0.005, 0.007, 0.01))),
                TeX("$f_b$"), TeX("$f_s$")),
       cex = 1,
       y.intersp = 1)



curve(fb_true, l, u, lwd = 2.5,
      col = 'brown', lty = 10, ylab = '')
curve(gb_test(x, fs_prop = 0), lwd = 2.5, add = TRUE,
      col = 'black', lty = 1)
curve(fs, lwd = 2.5, add = TRUE,
      col = 'skyblue', lty = 9)
legend('topright', col = c('brown', 'black', 'skyblue'),
       lty = c(10, 1, 9), bty = 'n', lwd =2.5,
       legend=c(TeX("$f_b$"), TeX("$q_b$"), TeX("$f_s$")),
       cex = 1,
       y.intersp = 1)



