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


###### SAFEGUARD ##################################################

file_bkg <- paste0('benchmark_toys/Cat',cat,'_',
                   scenario,'_mu0p0','.csv')
bkg_only_data <- read.csv('benchmark_toys/Cat0_HLHC_mu0p0.csv', header = FALSE)[,1] # change bkg-only data here

ni_bkg <- sapply(1:k, function(i){
  sum((bkg_only_data>bins[i])&(bkg_only_data<=bins[i+1]))
})
N_bkg <- sum(ni_bkg)

### DESIGNING qb ###############################################################

# qb likelihood:
qb_likelihood <- function(k, bin_ends, bin_counts)
{
  qb_mass <- integrate(function(x){
    x^k/((x-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  mu_i <- c()
  for(i in 1:(length(bin_ends)-1))
  {
    mu_i[i] <- N * integrate(f = function(x){
      (x^k/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bin_ends[i], upper = bin_ends[i+1])$value
  }
  return(sum(mu_i - log(mu_i)*bin_counts))
}

qb_pow_lik <- nlminb(start = 0.01,
                     objective = qb_likelihood,
                     lower = -Inf, upper = Inf,
                     bin_ends = bins, bin_counts = ni)$par

qb_mass_lik <- integrate(function(x){
  x^qb_pow_lik/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_lik <- function(x){
  (x^qb_pow_lik/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_lik
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
    x^k/((x-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  mu_i <- c()
  for(i in 1:(length(bins_left)-1))
  {
    mu_i <- c(mu_i,
              N * integrate(f = function(x){
                (x^k/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
              }, lower = bins_left[i], upper = bins_left[i+1])$value
    )}
  
  for(i in 1:(length(bins_right)-1))
  {
    mu_i <- c(mu_i, N * integrate(f = function(x){
      (x^k/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bins_right[i], upper = bins_right[i+1])$value
    )}
  
  return(sum(mu_i - log(mu_i)*c(ni_left, ni_right)))
}

qb_pow_sb <- nlminb(start = 0.01,
                    objective = qb_sb_likelihood,
                    lower = -Inf, upper = Inf,
                    bin_ends = bins, bin_counts = ni,
                    mid_region = c(120,130))$par

qb_mass_sb <- integrate(function(x){
  x^qb_pow_sb/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_sb <- function(x){
  (x^qb_pow_sb/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_sb
}


qb_pow_lik <- nlminb(start = 0.01,
                     objective = qb_likelihood,
                     lower = -Inf, upper = Inf,
                     bin_ends = bins, bin_counts = ni)$par

qb_mass_lik <- integrate(function(x){
  x^qb_pow_lik/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_lik <- function(x){
  (x^qb_pow_lik/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_lik
}

qb_pow_lik0 <- nlminb(start = 0.01,
                      objective = qb_likelihood,
                      lower = -Inf, upper = Inf,
                      bin_ends = bins, 
                      bin_counts = ni_bkg)$par

qb_mass_lik0 <- integrate(function(x){
  x^qb_pow_lik0/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_lik0 <- function(x){
  (x^qb_pow_lik0/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_lik0
}

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

curve(N*qb_lik0(x)*(u-l)/(k+1),
      col = alpha('blue', 0.5),
      lwd = 2.2, add = TRUE, lty = 4)

legend('topright',
       legend = TeX(c('BW','$\\hat{q}_b(MLE)$','$\\hat{q}_b(side-bands)$',
                      '$\\hat{q}_b(MLE -bkg)$')),
       bty = 'n',
       col = c('red', 'brown', 'orange', 'blue'),
       lwd = 2.2, lty = 1:4)