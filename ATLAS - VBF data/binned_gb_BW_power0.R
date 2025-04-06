rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

source('signal_search_functions.R')

mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110


dat <- read.csv('Data/VBF_Cat2.csv', header = FALSE) # change data here
ni <- dat[,2]
N <- sum(ni)
k <- nrow(dat)
bins <- seq(l, u, length.out = k+1)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

M_eps <- find_signal_region(Fs, 125, 1-1e-3)
M_lower <- M_eps[1]; M_upper <- M_eps[2]


# Defining q_b
bkg_loc <- 91.2
bkg_scale <- 2.49/2

obs <- generate_unbinned_data(ni, bins, seed = 1234)

BW_power_likelihood <- function(power)
{
  total_mass <- integrate(f = function(x){
    1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-power))
  }, l, u)$value
  
  fi <- sapply(obs, function(x) {
    (1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-power)))/total_mass
  })
  return(-sum(log(fi)))
}
(qb_pow <- nlminb(start = 0.01,
                  objective = BW_power_likelihood,
                  lower = 0.005, upper = 10)$par)

qb_ <- function(x) 1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-qb_pow))
qb_mass <- integrate(qb_, l, u)$value
# (pow_mass <- integrate(function(x) x^(-qb_pow), l, u)$value)
qb <- function(x) (1/((x-bkg_loc)^2 + bkg_scale^2) + (x^(-qb_pow)))/qb_mass

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 5)
res_sig_search <- binned_signal_search(lambda = lambda_seq,
                                       fs = fs, qb = qb,
                                       mu1 = mean1_in_gb, mu2 = mean2_in_gb, sd = sd_in_gb,
                                       search_region = c(l,u),
                                       count_region = c(120, 130),
                                       bin_mids = xi, bin_counts = ni)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c('lambda',
                              'eta_hat_binned', 'p-value_binned',
                              'S','significance')
knitr::kable(res_sig_search, 'pipe')

picture_l <- M_lower - 3; picture_u <- M_upper + 3
plot(xi, ni, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     col = alpha('black', alpha = 0.3),
     xlim = c(picture_l, picture_u))
mycols <- c('brown', 'skyblue' ,'red', 'orange', 'purple')
palette(mycols)
my_lty = 2:6

for(j in 1:length(lambda_seq))
{
  curve(N*((dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mean1_in_gb, sd = sd_in_gb) + 
              dtrunc(x, spec = 'norm', a = l, b = u,
                     mean = mean2_in_gb, sd = sd_in_gb))*lambda_seq[j] +
             (1-2*lambda_seq[j])*qb(x))*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 2.2,
        col = alpha(palette()[j], alpha = 0.4),
        lty = my_lty[j])
}

legend('bottomleft', col = mycols,
       lty = my_lty, bty = 'n', lwd =2.2,
       legend=c(TeX(sprintf(r'($g_b(\lambda = %f)$)', lambda_seq))),
       cex = 1,
       y.intersp = 1)
