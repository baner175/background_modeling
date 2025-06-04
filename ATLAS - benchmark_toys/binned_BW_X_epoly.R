rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)
library(knitr)
library(kableExtra)


source('signal_search_functions.R')

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110

# Defining BW parameters
bkg_loc <- 91.2
bkg_scale <- 2.49/2

### LOADING DATA ###############################################################
cat <- 0
scenario <- 'HLHC'
mu <- 1

mu_part_1 <- floor(mu); mu_part_2 <- (10*mu)%%10

file_name <- paste0('benchmark_toys/Cat',cat,'_',
                    scenario,'_mu',mu_part_1,'p',mu_part_2,'.csv')

obs <- read.csv(file_name, header = FALSE)[,1] # change data here

# binning the data
k <- 100 # number of bins
bins <- seq(l, u, length.out = k+1) # creating bin end-points
ni <- sapply(1:k, function(i){
  sum((obs>bins[i])&(obs<=bins[i+1]))
}) # creating bin counts
N <- sum(ni) # total number of events
xi <- (bins[1:k] + bins[2:(k+1)])/2 # creating bin mid-points

### DESIGNING qb ###############################################################

# qb likelihood:
qb_likelihood <- function(k, bin_ends, bin_counts)
{
  qb_mass <- integrate(function(x){
    exp(k*x)/((x-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  mu_i <- c()
  for(i in 1:(length(bin_ends)-1))
  {
    mu_i[i] <- N * integrate(f = function(x){
      (exp(k*x)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bin_ends[i], upper = bin_ends[i+1])$value
  }
  return(sum(mu_i - log(mu_i)*bin_counts))
}

qb_rate_lik <- nlminb(start = 0.01,
                      objective = qb_likelihood,
                      lower = -Inf, upper = Inf,
                      bin_ends = bins, bin_counts = ni)$par

qb_mass_lik <- integrate(function(x){
  exp(qb_rate_lik*x)/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_lik <- function(x){
  (exp(qb_rate_lik*x)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_lik
}


# qb sidebands likelihood
qb_sb_likelihood <- function(k, bin_ends, bin_counts,
                             mid_region)
{
  
  left_val <- mid_region[1]; right_val <- mid_region[2]
  
  bins_upper <- bin_ends[2:length(bin_ends)]
  bins_lower <- bin_ends[1:(length(bin_ends)-1)]
  
  bins_left <- bin_ends[bin_ends<=left_val]
  ni_left <- bin_counts[bins_upper<=left_val]
  
  bins_right <- bin_ends[bin_ends>=right_val]
  ni_right <- bin_counts[bins_lower>=right_val]
  
  N <- sum(ni_left) + sum(ni_right)
  qb_mass <- integrate(function(x){
    exp(k*x)/((x-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  mu_i <- c()
  for(i in 1:(length(bins_left)-1))
  {
    mu_i <- c(mu_i,
              N * integrate(f = function(x){
                (exp(k*x)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
              }, lower = bins_left[i], upper = bins_left[i+1])$value
    )}
  
  for(i in 1:(length(bins_right)-1))
  {
    mu_i <- c(mu_i, N * integrate(f = function(x){
      (exp(k*x)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bins_right[i], upper = bins_right[i+1])$value
    )}
  
  return(sum(mu_i - log(mu_i)*c(ni_left, ni_right)))
}

qb_rate_sb <- nlminb(start = 0.01,
                     objective = qb_sb_likelihood,
                     lower = -Inf, upper = Inf,
                     bin_ends = bins, bin_counts = ni,
                     mid_region = c(120,130))$par

qb_mass_sb <- integrate(function(x){
  exp(qb_rate_sb*x)/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_sb <- function(x){
  (exp(qb_rate_sb*x)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_sb
}

### SIGNAL SEARCH BEGINS #######################################################

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
M_eps <- find_signal_region(Fs, 125, 1-1e-3)
(M_lower <- M_eps[1]); (M_upper <- M_eps[2])

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 5)

# signal search using qb(MLE):
res_sig_search_qb_lik <- binned_signal_search(lambda = lambda_seq,
                                              fs = fs, qb = qb_lik,
                                              mu1 = mean1_in_gb, mu2 = mean2_in_gb,
                                              sd = sd_in_gb,
                                              search_region = c(l,u),
                                              bin_mids = xi, bin_counts = ni)
res_sig_search_qb_lik <- cbind(lambda_seq, t(res_sig_search_qb_lik))
res_sig_search_qb_lik <- data.frame(res_sig_search_qb_lik)


# signal search using qb(sidebands):
res_sig_search_qb_sb <- binned_signal_search(lambda = lambda_seq,
                                             fs = fs, qb = qb_sb,
                                             mu1 = mean1_in_gb, mu2 = mean2_in_gb,
                                             sd = sd_in_gb,
                                             search_region = c(l,u),
                                             bin_mids = xi, bin_counts = ni)
res_sig_search_qb_sb <- cbind(lambda_seq, t(res_sig_search_qb_sb))
res_sig_search_qb_sb <- data.frame(res_sig_search_qb_sb)

colnames(res_sig_search_qb_lik) <-
  colnames(res_sig_search_qb_sb) <- 
  c("$\\lambda$", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif', '$\\hat{S}$')

res_sig_search <- rbind(res_sig_search_qb_lik, 
                       res_sig_search_qb_sb)

caption <- paste0("Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2)
kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center') %>%
  column_spec(1:ncol(res_sig_search),
              extra_css = "padding-right: 30px;") %>%
  pack_rows("qb estimated with MLE", 1, nrow(res_sig_search_qb_lik)) %>%
  pack_rows("qb estimated with side bands", nrow(res_sig_search_qb_lik)+1,
            nrow(res_sig_search))

### MAKING PLOTS ###############################################################

# sensitivity analysis:
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
             (1-2*lambda_seq[j])*qb_lik(x))*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 2.2,
        col = alpha(palette()[j], alpha = 0.6),
        lty = my_lty[j])
}

legend('bottomleft', col = mycols,
       lty = 2:6, bty = 'n', lwd =2.2,
       legend=c(TeX(sprintf(r'($g_b(\lambda = %f)$)', lambda_seq))),
       cex = 1,
       y.intersp = 1)


# comparing q_b(MLE) and q_b(side-bands):
plot(x = xi, y = ni, pch = 16,
     col = alpha('black', 0.3),
     xlab = 'Energy (GeV)', ylab = 'Events')

curve(N*dtrunc(x, spec = 'cauchy', a = l, b = u,
               location = bkg_loc, scale = bkg_scale)*(u-l)/(k+1),
      col = alpha('red', 0.5),
      lwd = 2.2, add = TRUE, lty = 1)

curve(N*qb_lik(x)*(u-l)/(k+1),
      col = alpha('brown', 0.5),
      lwd = 2.2, add = TRUE, lty = 2)

curve(N*qb_sb(x)*(u-l)/(k+1),
      col = alpha('orange', 0.5),
      lwd = 2.2, add = TRUE, lty = 3)

legend('topright',
       legend = TeX(c('BW','$\\hat{q}_b(MLE)$','$\\hat{q}_b(side-bands)$')),
       bty = 'n',
       col = c('red', 'brown', 'orange'),
       lwd = 2.2, lty = 1:2)
