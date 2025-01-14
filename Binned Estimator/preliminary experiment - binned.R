rm(list = ls())
library("truncdist")
library(VGAM)

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
eta_true <- 0.03

n_breaks <- 100

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
bkg_samp <- rtrunc(n, a = l, b = u, spec = 'gamma',
                   rate = bkg_rate, shape = bkg_shape)

#generating data from signal
sig_samp <- rtrunc(n, spec = 'norm', a = l, b= u,
                   mean = mean_sig, sd = sd_sig)

u_mask <- runif(n)
samp <- ifelse(u_mask<eta_true, sig_samp, bkg_samp)

hs <- hist(samp, probability = TRUE, breaks = n_breaks)

binned_data <- hs$counts
binned_mids <- hs$mids

# desigining gb:

qb_beta <- 3.867429

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, spec = 'pareto', a = l, b = u,
                   scale = l, shape = qb_beta)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

search_signal <- function(lambda)
{
  norm_S <- integrate(function(t) {((fs(t)/gb_test(t, fs_prop = lambda) - 1)^2)*gb_test(t, fs_prop = lambda)}, l, u)$value |> sqrt()
  
  S1_vec <- sapply(binned_mids, function(x) {
    f_sig <- fs(x)
    g_b <- gb_test(x, fs_prop = lambda)
    return((f_sig/g_b -1)/norm_S)
  })
  theta <- sum(S1_vec*binned_data)/n
  # se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)
  
  # testing \eta = 0:
  eta_hat <- theta/norm_S
  
  # testing for signal:
  # (t_stat_theta <- theta/se_theta)
  # p_val <- pnorm(t_stat_theta, lower.tail = FALSE)
  
  return(eta_hat)
}

lambda_max <- 0.015
(lambda_seq <- seq(0, lambda_max, length.out = 4))
res_sig_search <- sapply(lambda_seq, search_signal)
res_sig_search <- cbind(lambda_seq, res_sig_search)
res_sig_search <- data.frame(res_sig_search)

mycols <- c('red', 'green', 'orange', 'purple')
palette(mycols)
my_lty = c(1,2,4,5)

picture_l <- 1.15 ; picture_u <- 1.4
plot(binned_mids, binned_data, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     xlim = c(picture_l, picture_u),
     main = paste0('# breaks = ',n_breaks))
arrows(x0 = binned_mids, y0 = binned_data + sqrt(binned_data),
       x1 = binned_mids, y1 = binned_data - sqrt(binned_data),
       angle = 90, length = 0.01, code = 3, col = 'red')
for(j in 1:length(lambda_seq))
{
  curve(n*(u-l)*gb_test(x, fs_prop = lambda_seq[j])/n_breaks, l, u,
        add = TRUE, lwd = 2.2,
        col = j,
        lty = my_lty[j])
}
abline(v = c(M_lower, M_upper), lwd = 2.2, col = 'black')
legend('topright', col = mycols,
       lty = my_lty, bty = 'n', lwd =2.2,
       legend=TeX(sprintf(r'($\lambda = %f (\hat{\eta}: %f)$)',
                            lambda_seq, round(as.numeric(res_sig_search[,2]), 4))),
       cex = 1.2,
       y.intersp = 1)
