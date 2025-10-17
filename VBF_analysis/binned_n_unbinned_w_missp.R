rm(list = ls())
library(truncdist)
library(VGAM)
library(latex2exp)
library(knitr)
library(kableExtra)
library(numDeriv)
library(nloptr)

# source('signal_search_functions.R')

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110
k <- 5e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

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

phys_data <- read.csv(file_name, header = FALSE)[,1]
N <- length(phys_data)
ni <- sapply(1:k, function(i){
  sum((phys_data>bin_ends[i])&(phys_data<=bin_ends[i+1]))
})

file_bkg <- paste0('benchmark_toys/Cat',cat,'_',
                   scenario,'_mu0p0','.csv')
bkg_data <- read.csv(file_bkg, header = FALSE)[,1]
M <- length(bkg_data)
mi <- sapply(1:k, function(i){
  sum((bkg_data>bin_ends[i])&(bkg_data<=bin_ends[i+1]))
})

spec = 'cauchy' # 'unif' or 'exp' or 'cauchy' 

# qb <- function(x) dtrunc(x, spec = spec, l, u, min = l, max = u)
# qb <- function(x) dtrunc(x, spec = spec, l, u, rate = 1/50)
qb <- function(x) dtrunc(x, spec = spec, l, u, location = 91.2, scale = 2.49/2)

norm_S <- integrate(function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- qb(x)
  S_val <- (fs/qb - 1)
  return((S_val^2)*qb)
}, l, u)$value |> sqrt()

S2 <- function(x) {
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- qb(x)
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S^2))
}

################################################################################
#-------------------------------------------------------------------------------
############ OUR BINNED METHOD WITH BKG DATA ###################################
S2_xi <- sapply(xi, S2)
theta0_hat_binned <- sum(ni*S2_xi)/N
delta0_hat_binned <- sum(mi*S2_xi)/M
eta_hat_binned <- (theta0_hat_binned - delta0_hat_binned)/(1-delta0_hat_binned)
cb_hat <- M/k

var_S2_F_hat_binned <- sum((S2_xi^2)*ni)/N - theta0_hat_binned^2
var_S2_Fb_hat_binned <- sum((S2_xi^2)*mi)/M - delta0_hat_binned^2

d_theta0_T_binned <- 1/(1-delta0_hat_binned)
d_delta0_T_binned <- (theta0_hat_binned-1)/((1-delta0_hat_binned)^2)


test_num_binned <- sqrt(M*N)*(eta_hat_binned - 0)

denom1_binned <- (M/((1-delta0_hat_binned)^2)) * var_S2_F_hat_binned

denom2_binned <- N*(((theta0_hat_binned-1)^2)/((1-delta0_hat_binned)^4)) * var_S2_Fb_hat_binned


test_denom_binned <- sqrt(denom1_binned + denom2_binned)

eta_test_stat_binned <- test_num_binned/test_denom_binned
p_val_binned <- pnorm(eta_test_stat_binned,lower.tail = FALSE)

unbiased_row_binned <- data.frame(
  'unbiased test - binned',
  eta_hat_binned, 
  p_val_binned, 
  max(qnorm(p_val_binned, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_hat_binned*N, 2), 0))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
############ OUR UNBINNED METHOD WITH BKG DATA ###################################
S2_phys_vec <- sapply(phys_data, S2)
S2_bkg_vec <- sapply(bkg_data, S2)

theta0_hat_unbinned <- mean(S2_phys_vec)
delta0_hat_unbinned <- mean(S2_bkg_vec)
eta_hat_unbinned <- (theta0_hat_unbinned - delta0_hat_unbinned)/(1-delta0_hat_unbinned)

var_S2_F_hat_unbinned <- mean(S2_phys_vec^2) - theta0_hat_unbinned^2
var_S2_Fb_hat_unbinned <- mean(S2_bkg_vec^2) - delta0_hat_unbinned^2

d_theta0_T_unbinned <- 1/(1-delta0_hat_unbinned)
d_delta0_T_unbinned <- (theta0_hat_unbinned-1)/((1-delta0_hat_unbinned)^2)


test_num_unbinned <- sqrt(M*N)*(eta_hat_unbinned - 0)

denom1_unbinned <- (M/((1-delta0_hat_unbinned)^2)) * var_S2_F_hat_unbinned

denom2_unbinned <- N*(((theta0_hat_unbinned-1)^2)/((1-delta0_hat_unbinned)^4)) * var_S2_Fb_hat_unbinned


test_denom_unbinned <- sqrt(denom1_unbinned + denom2_unbinned)

eta_test_stat_unbinned <- test_num_unbinned/test_denom_unbinned
p_val_unbinned <- pnorm(eta_test_stat_unbinned,lower.tail = FALSE)

unbiased_row_unbinned <- data.frame(
  'unbiased test - unbinned',
  eta_hat_unbinned, 
  p_val_unbinned, 
  max(qnorm(p_val_unbinned, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_hat_unbinned*N, 2), 0))

#-------------------------------------------------------------------------------
################################################################################

######## Creating results table ################################################

colnames(unbiased_row_binned) <- colnames(unbiased_row_unbinned) <- 
  c("Test", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif',
    '$\\hat{S}$')

res_sig_search <- rbind(unbiased_row_binned,
                        unbiased_row_unbinned)

caption <- paste0("Binned (", k," bins) and Unbinned Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2, '<br> Bkg - ', spec)

kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center')

################################################################################
######### Plotting Densities ###################################################

k_graph <- 100
bin_ends_graph <- seq(l, u, length.out = k_graph + 1)
xi_graph <- (bin_ends_graph[1:k_graph] + bin_ends_graph[2:(k_graph+1)])/2
ni_graph <- sapply(1:k_graph, function(i){
  sum((phys_data>bin_ends_graph[i])&(phys_data<=bin_ends_graph[i+1]))
})
mi_graph <- sapply(1:k_graph, function(i){
  sum((bkg_data>bin_ends_graph[i])&(bkg_data<=bin_ends_graph[i+1]))
})


plot(x = xi_graph, y = mi_graph, pch = 16, col = 'grey',
     xlab = 'Energy(GeV)', ylab = 'Bkg counts',
     main = paste0('Misspecified density - ', spec,' on bkg'))

curve(M*qb(x)*(u-l)/(k_graph+1), l, u, add = TRUE,
      col = ggplot2::alpha('black', 0.6), lwd = 2, lty = 2)
