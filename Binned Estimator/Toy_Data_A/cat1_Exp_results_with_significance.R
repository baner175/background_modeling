rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

mean_sig <- 125; sd_sig <- 2
eps <- 1e-3
u <- 160; l <- 110

dat <- read.csv('Data/histogram_data_cat1.csv', header = TRUE)

ni <- dat$Content
N <- sum(ni)
k <- nrow(dat)
bins <- seq(l, u, 0.1)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# Generating pseudo unbinned data to generate KDE
set.seed(seed = 12345)
obs_unbinned <- c()
for(i in 1:k)
{
  obs_unbinned <- c(obs_unbinned, runif(ni[i], bins[i], bins[i+1]))
}

# obs_normalised <- scale(obs_unbinned)
# kde_ <- kdensity(obs_normalised, bw = 0.1)
# kde_unnormed<- function(t) kde_((t-mean(obs_unbinned))/sd(obs_unbinned))
# kde_tot <- integrate(kde_unnormed, l, u)$value
# kde <- function(x) kde_unnormed(x)/kde_tot
# plot(x = xi, y = ni, pch = 16)
# curve(N*kde(x)*(u-l)/(k+1), col = 'blue', add = TRUE, lwd = 2)

########################################################################

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

# Defining g_b
## likelihood using pseudo-unbinned data
qb_likelihood <- function(pars)
{
  rate <- pars
  
  fi <- sapply(obs_unbinned, function(t)
  {
    dtrunc(t, rate = rate,
           spec = 'exp',
           a = l, b = u)
  })
  return(-sum(log(fi)))
}

# (res <- nlminb(start = 0.01,
#                objective = qb_likelihood,
#                lower = 0.005, upper = 0.3)$par) # TAKES TOO LONG TO RUN SO USING THE VALUE FROM CAT0
qb_rate <- 0.052

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, a = l, b= u, spec = 'exp', 
                   rate = qb_rate)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}


(N_mid <- sum(ni[xi<=130 & xi>=120]))

search_signal <- function(lambda)
{
  norm_S <- integrate(function(t) {((fs(t)/gb_test(t, fs_prop = lambda) - 1)^2)*gb_test(t, fs_prop = lambda)}, l, u)$value |> sqrt()
  
  # Testing using the binned data:
  S1_vec <- sapply(xi, function(x) {
    f_sig <- fs(x)
    g_b <- gb_test(x, fs_prop = lambda)
    return((f_sig/g_b -1)/norm_S)
  })
  theta_hat_binned <- sum(S1_vec*ni)/N
  
  # estimating eta:
  eta_hat_binned <- theta_hat_binned/norm_S
  
  # testing for signal:
  theta_check <- mean(S1_vec*ni)
  t_stat_theta <- k*theta_check/sqrt(sum(ni*S1_vec^2))
  p_val_binned <- pnorm(t_stat_theta, lower.tail = FALSE)
  S_hat <- N_mid*eta_hat_binned
  B_hat <- N_mid*(1-eta_hat_binned)
  signif <- S_hat/sqrt(B_hat)
  
  return(c(eta_hat_binned, p_val_binned, 
           round(S_hat, 2), round(signif, 3)))
  
}


lambda_max <- 0.005
(lambda_seq <- seq(0, lambda_max, length.out = 4))
res_sig_search <- sapply(lambda_seq, search_signal)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c('lambda',
                              'eta_hat_binned', 'p-value_binned',
                              'S','significance')
res_sig_search


picture_l <- M_lower - 3; picture_u <- M_upper + 3
plot(xi, ni, pch = 16,
     col = alpha('black', alpha = 0.3),
     ylab="Events",xlab="Energy (Gev)",
     xlim = c(picture_l, picture_u))

# arrows(x0 = xi, y0 = ni + sqrt(ni),
#        x1 = xi, y1 = ni - sqrt(ni),
#        angle = 90, length = 0.01, code = 3, col = 'red')

abline(v = c(M_lower, M_upper), col = 'black', lwd = 2.2)

mycols <- c('red', 'green', 'orange', 'purple')
palette(mycols)
my_lty = c(3,4,5,6)

for(j in 1:length(lambda_seq))
{
  curve(N*gb_test(x, fs_prop = lambda_seq[j])*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 2.2,
        col = j,
        lty = my_lty[j])
}

legend('bottomleft', col = mycols,
       lty = my_lty, bty = 'n', lwd =2.2,
       legend=c(TeX(sprintf(r'($\lambda = %f (\hat{\eta}_{(b)}: %f , p.val_{(b)} : %f)$)', lambda_seq,
                            round(as.numeric(res_sig_search[,2]), 4),
                            round(as.numeric(res_sig_search[,3]), 4)))),
       cex = 1,
       y.intersp = 1)
title(main = 'Toy Data A - Cat 1')

knitr::kable(res_sig_search[-1,], 'pipe')
