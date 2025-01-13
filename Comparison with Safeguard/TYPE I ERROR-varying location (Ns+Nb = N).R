rm(list = ls())

library("truncdist")
library("VGAM")

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
curve(fb_true, lwd = 2, col = 'red', add = TRUE, lty = 2)
curve(dtrunc(x, spec = 'pareto', a = l, b = u,
             shape = qb_beta, scale = l),
      lwd = 2, col = 'brown', add = TRUE)

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
curve(gb_sf, col = 'cyan', lwd = 2, from = l, to = u, add = TRUE)
curve(fs, col = 'blue', lwd = 2, from = l, to = u, add = TRUE)
curve(fb_true, col = 'brown', lwd = 2, from = l, to = u, add = TRUE)
legend('topright', col = c('cyan', 'brown', 'blue'),
       legend = c('proposed background from safeguard', 'true background', 'signal'),
       bty = 'n', lty = 1, lwd = 2)
abline(v = c(M_lower, M_upper), lwd = 2, col = 'black')

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

curve(gb_test(x, fs_prop = 0), col = 'blue', lty = 1, lwd = 2,
      from = l, to = u, ylab = 'gb_test')

curve(gb_test(x, fs_prop = 0.002), col = 'blue', lty = 2, lwd = 2,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.005), col = 'blue', lty = 2, lwd = 2,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.007), col = 'blue', lty = 2, lwd = 2,
      add = TRUE)
curve(gb_test(x, fs_prop = 0.01), col = 'blue', lty = 2, lwd = 2,
      add = TRUE)
curve(gb_sf, col = 'cyan', lwd = 2, lty = 1, add = TRUE)

(norm_S_1 <- integrate(function(t) 
{((fs(t)/gb_test(t, fs_prop = 0.002) - 1)^2)*gb_test(t, fs_prop = 0.002)},
l, u)$value |> sqrt())

(norm_S_2 <- integrate(function(t) 
{((fs(t)/gb_test(t, fs_prop = 0.005) - 1)^2)*gb_test(t, fs_prop = 0.005)},
l, u)$value |> sqrt())

(norm_S_3 <- integrate(function(t) 
{((fs(t)/gb_test(t, fs_prop = 0.007) - 1)^2)*gb_test(t, fs_prop = 0.007)},
l, u)$value |> sqrt())

(norm_S_4 <- integrate(function(t) 
{((fs(t)/gb_test(t, fs_prop = 0.01) - 1)^2)*gb_test(t, fs_prop = 0.01)},
l, u)$value |> sqrt())


###################################################################
######################## SIMULATION STARTS ########################
###################################################################

set.seed(12344)
B <- 1e5
p_vals <- c()
test_stat <- c()
eta_hat <- c()

time_now <- Sys.time()
for(j in 1:B)
{
  n_sim <- rpois(1, lambda = n)
  # n_sim <- n
  samp_sim <- rtrunc(n_sim,a=l,b=u,spec="gamma",
                     rate=bkg_rate,shape=bkg_shape)
  
  ######################### SAFEGUARD #########################
  
  sf_sim <- nlminb(start = 0.01,
                   objective = mix_mod_sf,
                   lower = 0,
                   upper = 1,
                   data = samp_sim)
  
  # adding a minus sign since mix_mod returns negative log likelihood
  ll0 <- -mix_mod_sf(pars = 0, data = samp_sim)
  ll1 <- -sf_sim$objective
  test_stat_sf <- -2*(ll0-ll1)
  p_val_sf <- 0.5*pchisq(test_stat_sf, df = 1, lower.tail = FALSE)
  eta_hat_sf <- sf_sim$par
  
  ######################## OUR TEST ############################
  
  ####################### lambda = 0.002 ########################
  
  S1_vec <- (fs(samp_sim)/gb_test(samp_sim, fs_prop = 0.002) - 1)/norm_S_1
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n_sim)
  
  t_stat_1 <- theta/se_theta
  p_val_1 <- pnorm(t_stat_1, lower.tail = FALSE)
  eta_hat_1 <- theta/norm_S_1
  
  
  ####################### lambda = 0.005 ########################
  
  S1_vec <- (fs(samp_sim)/gb_test(samp_sim, fs_prop = 0.005) - 1)/norm_S_2
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n_sim)
  
  t_stat_2 <- theta/se_theta
  p_val_2 <- pnorm(t_stat_2, lower.tail = FALSE)
  eta_hat_2 <- theta/norm_S_2
  
  ####################### lambda = 0.007 ########################
  
  S1_vec <- (fs(samp_sim)/gb_test(samp_sim, fs_prop = 0.007) - 1)/norm_S_3
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n_sim)
  
  t_stat_3 <- theta/se_theta
  p_val_3 <- pnorm(t_stat_3, lower.tail = FALSE)
  eta_hat_3 <- theta/norm_S_3
  
  ####################### lambda = 0.01 ########################
  
  S1_vec <- (fs(samp_sim)/gb_test(samp_sim, fs_prop = 0.01) - 1)/norm_S_4
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n_sim)
  
  t_stat_4 <- theta/se_theta
  p_val_4 <- pnorm(t_stat_4, lower.tail = FALSE)
  eta_hat_4 <- theta/norm_S_4
  
  ###############################################################
  
  p_vals <- rbind(p_vals, 
                  c(p_val_sf, p_val_1, p_val_2, p_val_3, p_val_4))
  test_stat <- rbind(test_stat, 
                     c(test_stat_sf, t_stat_1, t_stat_2, t_stat_3, t_stat_4))
  eta_hat <- rbind(eta_hat, 
                   c(eta_hat_sf, eta_hat_1, eta_hat_2, eta_hat_3, eta_hat_4))
  
  
  cat(sprintf('\rIteration: %d/%d',j,B))
}
Sys.time() - time_now # Took about 14 hours

p_vals <- data.frame(p_vals)
eta_hat <- data.frame(eta_hat)
test_stat <- data.frame(test_stat)

colnames(p_vals) <- colnames(eta_hat) <- colnames(test_stat) <- 
  c('Safeguard', 'fs_prop=0.002', 'fs_prop=0.005', 'fs_prop=0.007', 'fs_prop=0.01')

# write.csv(p_vals,
#           file = 'Simulation Results/Type 1 error (Ns+Nb=N) p-values.csv',
#           row.names = FALSE)
# 
# write.csv(eta_hat,
#           file = 'Simulation Results/Type 1 error (Ns+Nb=N) eta_hat.csv',
#           row.names = FALSE)
# 
# write.csv(test_stat,
#           file = 'Simulation Results/Type 1 error (Ns+Nb=N) test_stat.csv',
#           row.names = FALSE)


p_vals <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) p-values.csv', header = TRUE)
eta_hat <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) eta_hat.csv', header = TRUE)
test_stat <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) test_stat.csv', header = TRUE)

# Type 1 error for safeguard:
mean(p_vals$Safeguard<=0.05) # 0.2417 # corresponding delta: 0.01340076

# checking delta:
sapply(c(0.002, 0.005, 0.007, 0.01), function(lam){
  normS <- integrate(function(t) {((fs(t)/gb_test(t, fs_prop = lam) - 1)^2)*gb_test(t, fs_prop = lam)},
                     l, u)$value |> sqrt()
  
  del <- integrate(function(t) (fs(t)/gb_test(t, fs_prop = lam) - 1)*fb_true(t),
                   l, u)$value/normS
  return(del)
}
)

# Type 1 error for fs_prop = 0.002:
mean(p_vals$fs_prop.0.002<=0.05) # 0.29077 # corresponding delta: 0.016249060

# Type 1 error for fs_prop = 0.005:
mean(p_vals$fs_prop.0.005<=0.05) # 0.12836 # corresponding delta: 0.007890316

# Type 1 error for fs_prop = 0.007:
mean(p_vals$fs_prop.0.007<=0.05) # 0.06422 # corresponding delta: 0.002424086

# Type 1 error for fs_prop = 0.01:
mean(p_vals$fs_prop.0.01<=0.05) # 0.01748 # corresponding delta: -0.005624272
