rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

mean_sig <- 125; sd_sig <- 2
eps <- 1e-3
u <- 160; l <- 110

dat <- read.csv('Data/histogram_data_cat3.csv', header = TRUE)

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

obs_normalised <- scale(obs_unbinned)
kde_ <- kdensity(obs_normalised, bw = 0.1)
kde_unnormed<- function(t) kde_((t-mean(obs_unbinned))/sd(obs_unbinned))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
plot(x = xi, y = ni, pch = 16)
curve(N*kde(x)*(u-l)/(k+1), col = 'blue', add = TRUE, lwd = 2)

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

# Defining q_b and g_b

bkg_loc <- 91.2
bkg_scale <- 2.49/2

qb <- function(x)
{
  dtrunc(x, a = l, b= u, spec = 'cauchy', 
         location = bkg_loc, 
         scale = bkg_scale)
}

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, a = l, b= u, spec = 'cauchy', 
                   location = bkg_loc, 
                   scale = bkg_scale)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}


# Area under gb
integrate(gb_test,l,u)

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
  
  
  # Testing using the pseudo-ubinned data:
  S1_vec <- sapply(obs_unbinned, function(x) {
    f_sig <- fs(x)
    g_b <- gb_test(x, fs_prop = lambda)
    return((f_sig/g_b -1)/norm_S)
  })
  theta_hat_unbinned <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta_hat_unbinned^2)/N)
  
  # testing \eta = 0:
  eta_hat_unbinned <- theta_hat_unbinned/norm_S
  
  # testing for signal:
  t_stat_theta <- theta_hat_unbinned/se_theta
  p_val_unbinned <- pnorm(t_stat_theta, lower.tail = FALSE)
  
  return(c(eta_hat_binned, p_val_binned,
           eta_hat_unbinned, p_val_unbinned))
  
}


lambda_max <- 0.005
(lambda_seq <- seq(0, lambda_max, length.out = 4))
res_sig_search <- sapply(lambda_seq, search_signal)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c('lambda',
                              'eta_hat_binned', 'p-value_binned',
                              'eta_hat_unbinned', 'p-value_unbinned')
res_sig_search


picture_l <- M_lower - 3; picture_u <- M_upper + 3
plot(xi, ni, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     type = 'n',
     xlim = c(picture_l, picture_u))
# arrows(x0 = xi, y0 = ni + sqrt(ni),
#        x1 = xi, y1 = ni - sqrt(ni),
#        angle = 90, length = 0.01, code = 3, col = 'red')

curve(N*kde(x)*(u-l)/(k+1), lwd = 2.2, col = 'black', add = TRUE, lty = 1)
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

legend('topright', col = c('black', mycols),
       lty = c(1, 2, my_lty), bty = 'n', lwd =2.2,
       legend=c('KDE',
                TeX(sprintf(r'($\lambda = %f (\hat{\eta}_{(b)}: %f , p.val_{(b)} : %f, \hat{\eta}_{(ub)}: %f , p.val_{(ub)} : %f)$)', lambda_seq,
                            round(as.numeric(res_sig_search[,2]), 4),
                            round(as.numeric(res_sig_search[,3]), 4),
                            round(as.numeric(res_sig_search[,4]), 4),
                            round(as.numeric(res_sig_search[,5]), 4)))),
       cex = 1,
       y.intersp = 1)
title(main = 'Toy Data C - Cat 3')


