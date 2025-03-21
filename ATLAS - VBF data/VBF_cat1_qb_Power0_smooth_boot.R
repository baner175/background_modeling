rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)
library(Hmisc)

source('necessary_functions.R')

mean_sig <- 125; sd_sig <- 3
eps <- 1e-3
u <- 160; l <- 110
dat <- read.csv('Data/VBF_Cat1.csv', header = FALSE)

ni <- dat[,2]
N <- sum(ni)
k <- nrow(dat)
bins <- seq(l, u, length.out = k+1)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# generating pseudo unbinned data
set.seed(123456)
obs <- c()
for(i in 1:(length(bins))-1)
{
  obs <- c(obs, runif(ni[i], bins[i], bins[i+1]))
}

# Defining signal density and calculating signal region

fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

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

# likelihood using pseudo-unbinned data
qb_likelihood <- function(pars)
{
  shape <- pars
  
  fi <- sapply(obs, function(t)
  {
    dtrunc(t, shape = shape, scale = l,
           spec = 'pareto',
           a = l, b = u)
  })
  return(-sum(log(fi)))
}

(res <- nlminb(start = 0.01,
               objective = qb_likelihood,
               lower = 0.005, upper = 10)$par)

qb_shape <- 3.35

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, a = l, b= u, spec = 'pareto', 
                   shape = qb_shape, scale = l) # change qb here
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

plot(xi, ni, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     col = alpha('black', alpha = 0.3))

# arrows(x0 = xi, y0 = ni + sqrt(ni),
#        x1 = xi, y1 = ni - sqrt(ni),
#        angle = 90, length = 0.01, code = 3,
#        col = alpha('red', alpha = 0.3))

curve(N*gb_test(x, fs_prop = 0)*(u-l)/(k+1),
      l, u, add = TRUE, lwd = 2.2,
      col = 'red',
      lty = 1)

curve(N*fs(x)*(u-l)/(k+1),
      l, u, add = TRUE, lwd = 2.2,
      col = 'skyblue',
      lty = 2)
abline(v = c(M_lower, M_upper), col = 'black', lwd = 2)

(N_mid <- sum(ni[xi<=130 & xi>=120]))
n_basis <- 2

search_signal <- function(lambda, boot = 1e5, seed = 12345)
{
  gb <- function(x) gb_test(x, fs_prop = lambda)
  Gb <- function(x, ...)
  {
    integrate(function(t) gb(t, ...), l, x)$value
  }
  
  Gb <- Vectorize(Gb)
  norm_S <- integrate(function(t) {((fs(t)/gb(t) - 1)^2)*gb(t)}, l, u)$value |> sqrt()
  
  basis <- construct_basis(n_basis = n_basis,
                           sig_density = fs,
                           proposal_bkg = gb,
                           limits = c(l,u))
  S1 <- basis[[2]]
  T_basis <- basis[-c(1,2)]
  tau <- sapply(T_basis, function(f) sapply(obs, f) |> mean())
  fm_null <- function(x)
  {
    T_vec <- sapply(T_basis, function(f) f(x))
    return(gb(x)*(1 + crossprod(tau, T_vec)[1,1]))
  }  
  fm_null <- Vectorize(fm_null)
  t_stat_boot <- c()
  
  ## calculating test statistic using 
  S1_vec <- sapply(xi, function(x) {
    f_sig <- fs(x)
    g_b <- gb_test(x, fs_prop = lambda)
    return((f_sig/g_b -1)/norm_S)
  })
  
  # bootstrapping
  set.seed(seed)
  
  u_obs <- Gb(obs)
  extrap <- approxExtrap(obs, fm_null(obs)/gb(obs), xout = c(l, obs, u))
  dG_GF <- approxfun(extrap$x,extrap$y, rule=2) # d(Gb(x),; Gb, Fm_hat)
  extrap <- approxExtrap(u_obs, fm_null(obs)/gb(obs), xout = c(0, u_obs, 1))
  du_GF <- approxfun(extrap$x,extrap$y, rule=2)
  M <- max(sapply(u_obs, du_GF))
  n_sim <- round(1.5*N*M)
  fs_prop <- lambda
  
  for(b in 1:boot)
  {
    ## sampling from Gb
    xFs1 <- rtrunc(n_sim, spec = 'norm', a = l, b = u,
                   mean = mean1_in_gb, sd = sd_in_gb)
    xFs2 <- rtrunc(n_sim, spec = 'norm', a = l, b = u,
                   mean = mean2_in_gb, sd = sd_in_gb)
    vG_sim1 <- runif(n_sim)
    xG_bump <- ifelse(vG_sim1<0.5, xFs1, xFs2)
    xQ <- rtrunc(n_sim, a = l, b= u, spec = 'pareto', 
                 shape = qb_shape, scale = l) # change qb here
    vG_sim2 <- runif(n_sim)
    xG <- ifelse(vG_sim2<1-2*fs_prop, xQ, xG_bump) # sample from Gb
    
    ## sampling from fm_null
    dG_GF.xG <- sapply(xG, dG_GF)
    v_sim <- runif(n_sim)
    xF <- xG[v_sim*M<dG_GF.xG]
    xF <- xF[1:N]
    ni_boot <- hist(xF, breaks = bins, plot = FALSE)$count
    
    theta_check_boot <- mean(S1_vec*ni_boot)
    
    t_stat_boot[b] <- k*theta_check_boot/sqrt(sum(ni_boot*S1_vec^2))
    cat((sprintf('\rLambda: %f; Iteration: %d/%d', lambda, b, boot)))
  }
  
  
  theta_hat_binned <- sum(S1_vec*ni)/N
  
  # estimating eta:
  eta_hat_binned <- theta_hat_binned/norm_S
  
  # testing for signal:
  theta_check <- mean(S1_vec*ni)
  t_stat_theta <- k*theta_check/sqrt(sum(ni*S1_vec^2))
  
  p_val_boot <- mean(t_stat_boot>t_stat_theta)
  S_hat <- N_mid*eta_hat_binned
  B_hat <- N_mid*(1-eta_hat_binned)
  signif <- S_hat/sqrt(B_hat)
  return(c(eta_hat_binned, p_val_boot, round(S_hat, 2), round(signif, 3)))
}

lambda_max <- 0.005
(lambda_seq <- seq(0, lambda_max, length.out = 4))
res_sig_search <- sapply(lambda_seq, search_signal)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c('lambda',
                              'eta_hat_binned', 'smooth booted p-value',
                              'S','significance')
res_sig_search

options(digits = 5)

picture_l <- M_lower - 3; picture_u <- M_upper + 3
plot(xi, ni, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     col = alpha('black', alpha = 0.3),
     xlim = c(picture_l, picture_u))

# arrows(x0 = xi, y0 = ni + sqrt(ni),
#        x1 = xi, y1 = ni - sqrt(ni),
#        angle = 90, length = 0.01, code = 3,
#        col = alpha('red', alpha = 0.3))


mycols <- c('red', 'green', 'orange', 'purple')
palette(mycols)
my_lty = c(3,4,5,6)

for(j in 1:length(lambda_seq))
{
  curve(N*gb_test(x, fs_prop = lambda_seq[j])*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 2.2,
        col = alpha(palette()[j], alpha = 0.4),
        lty = my_lty[j])
}

legend('bottomleft', col =  mycols,
       lty = my_lty, bty = 'n', lwd =2.2,
       legend=c(TeX(sprintf(r'($\lambda = %f (\hat{\eta}_{(b)}: %f)$)', lambda_seq,
                            round(as.numeric(res_sig_search[,2]), 4)))),
       cex = 1,
       y.intersp = 1)

knitr::kable(res_sig_search, 'pipe')
