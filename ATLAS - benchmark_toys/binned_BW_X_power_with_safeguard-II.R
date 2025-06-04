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
k <- 100
bins <- seq(l, u, length.out = k+1)
ni <- sapply(1:k, function(i){
  sum((obs>bins[i])&(obs<=bins[i+1]))
})
N <- sum(ni)
xi <- (bins[1:k] + bins[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

###### SAFEGUARD ##################################################

# Here we shall estimate beta in qb for only mu = 125, and then we'll estimate
# epsilon for varying mu between 120 to 130 by fixng beta at the 
# aforementioned estimate

file_bkg <- paste0('benchmark_toys/Cat',cat,'_',
                   scenario,'_mu0p0','.csv')
bkg_only_data <- read.csv(file_bkg, header = FALSE)[,1] # change bkg-only data here

ni_bkg <- sapply(1:k, function(i){
  sum((bkg_only_data>bins[i])&(bkg_only_data<=bins[i+1]))
})
N_bkg <- sum(ni_bkg)

### DESIGNING qb ###############################################################

# WITH THE BKG ONLY DATA, FOR THE TIME BEING 
# WE ARE USING AN UNBINNED LIKELIHOOD

sf_beta_eps_likelihood <- function(par, mu, bin_ends,  bin_counts)
{
  eps <- par[1]; power <- par[2]
  
  qb_mass <- integrate(function(t){
    (t^power)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  N_total <- sum(bin_counts)
  bi <- bin_ends
  mu_i <- c()
  for(i in 1:(length(bi)-1))
  {
    mu_i[i] <- N_total * integrate(f = function(x){
      eps*dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu, sd = sd_sig) + 
        (1-eps)*((x^power)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bi[i], upper = bi[i+1])$value
  }
  
  cat(sprintf("\rEvaluating at mu: %f", mu))
  return(sum(mu_i - log(mu_i)*bin_counts))
}

qb_pow <-  nlminb(start = c(0.01, 0.01),
                  objective = sf_beta_eps_likelihood,
                  lower = c(0, -Inf), upper = c(1, Inf),
                  mu = mean_sig, bin_ends = bins,
                  bin_counts = ni_bkg)$par[2]

sf_eps_likelihood <- function(par, mu, bin_ends,  bin_counts)
{
  eps <- par; power <- qb_pow
  
  qb_mass <- integrate(function(t){
    (t^power)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  N_total <- sum(bin_counts)
  bi <- bin_ends
  mu_i <- c()
  for(i in 1:(length(bi)-1))
  {
    mu_i[i] <- N_total * integrate(f = function(x){
      eps*dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mu, sd = sd_sig) + 
        (1-eps)*((x^power)/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, lower = bi[i], upper = bi[i+1])$value
  }
  
  cat(sprintf("\rEvaluating at mu: %f", mu))
  return(sum(mu_i - log(mu_i)*bin_counts))
}

mu_seq <- seq(120, 130, length.out = 21)
sf_res <- sapply(mu_seq, function(t)
{
  nlminb(start = 0.01,
         objective = sf_eps_likelihood,
         lower = 0, upper = 1,
         mu = t, bin_ends = bins,
         bin_counts = ni_bkg)$par
})

(gb_sf_eps <- max(sf_res))

# qb used for safeguard
qb_ <- function(x) (x^(qb_pow))/((x-bkg_loc)^2 + bkg_scale^2)
qb_mass <- integrate(qb_, l, u)$value
qb <- function(x) ((x^(qb_pow))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass

gb_sf <- function(x){
  gb_sf_eps*fs(x) + (1-gb_sf_eps)*qb(x)
}

norm_S_gb_sf <- integrate(function(t) ((fs(t)/gb_sf(t) - 1)^2)*gb_sf(t),
                          l, u)$value |> sqrt()

sf_phys_likelihood <- function(eta, bin_ends, bin_counts)
{
  N_total <- sum(bin_counts)
  bi <- bin_ends
  mu_i <- c()
  for(i in 1:(length(bi)-1))
  {
    mu_i[i] <- N_total * integrate(f = function(x){
      eta*dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig) + 
        (1-eta)*(gb_sf_eps*dtrunc(x, spec = 'norm', a = l, b = u,
                                  mean = mean_sig, sd = sd_sig) + 
                   (1-gb_sf_eps)*(((x^(qb_pow))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass))
    }, lower = bi[i], upper = bi[i+1])$value
  }
  return(sum(mu_i - log(mu_i)*bin_counts))
}

safeguard_res <- nlminb(start = 0,
                        objective = sf_phys_likelihood,
                        lower = 0, upper = 1,
                        bin_ends = bins, bin_counts = ni)

# estimated proportion from safeguard:
(eta_sf <- safeguard_res$par)

ll0 <- -sf_phys_likelihood(eta = 0,bin_ends = bins,
                           bin_counts = ni)
ll1 <- -safeguard_res$objective

# safeguard p-value:
(p_val_sf <- 0.5*pchisq(-2*(ll0-ll1),
                        df = 1,
                        lower.tail = FALSE))

sf_row <- data.frame(
  'safeguard', eta_sf, p_val_sf, qnorm(p_val_sf, lower.tail = FALSE) |> round(1),
  round(eta_sf*N, 2))

############ OUR METHOD #########################################

M_eps <- find_signal_region(Fs, mean_sig, 1-1e-3)
(M_lower <- M_eps[1]); (M_upper <- M_eps[2])

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 5)
res_sig_search <- binned_signal_search(lambda = lambda_seq,
                                       fs = fs, qb = qb,
                                       mu1 = mean1_in_gb, mu2 = mean2_in_gb, sd = sd_in_gb,
                                       search_region = c(l,u),
                                       bin_mids = xi, bin_counts = ni)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
res_sig_search[,1] <- as.character(res_sig_search[,1])

######### Using bkg-only data to estimate delta ################

S1_vec_sf <- sapply(bkg_only_data, function(t){
  (fs(t)/gb_sf(t) - 1)/norm_S_gb_sf
})
delta_hat_sf <- mean(S1_vec_sf)

delta_se_sf <- sqrt((mean(S1_vec_sf^2) - delta_hat_sf^2)/N_bkg)

lambda_i <- seq(0, 0.1, length.out = 1e2)

delta_res <- sapply(lambda_i, function(t){
  unbinned_delta_estimation(lambda = t, fs = fs, qb = qb,
                            mu1 = mean1_in_gb, mu2 = mean2_in_gb,
                            sd = sd_in_gb, search_region = c(l,u),
                            bkg_data = bkg_only_data)
})

delta_i <- delta_res[1,]
delta_i_se <- delta_res[2,]

plot(x = lambda_i, y = delta_i, type = 'l', lwd = 2,
     col = 'blue', ylim = range(delta_i+c(-3,3)*delta_i_se*1.96),
     xlab = TeX('$\\lambda$'),
     ylab = TeX('$\\hat{\\delta}$'))

polygon(c(lambda_i, rev(lambda_i)), 
        c(delta_i - delta_i_se*1.96, rev(delta_i + delta_i_se*1.96)), 
        col = rgb(0.1, 0.4, 0.8, 0.2), border = NA)
abline(h = 0, col = 'red', lwd = 2, lty = 2)

delta_fun <- approxfun(x = lambda_i, y = delta_i, rule = 2)
lambda0 <- uniroot(delta_fun, c(0, 0.1))$root
abline(v = lambda0, col = 'red', lty = 2, lwd = 2)
points(x = lambda0, y = 0, pch = 16, col = 'red')

unb_res <- binned_signal_search_1(lambda = lambda0, 
                                  fs = fs, qb = qb,
                                  mu1 = mean1_in_gb, 
                                  mu2 = mean2_in_gb, sd = sd_in_gb,
                                  search_region = c(l,u),
                                  bin_mids = xi, bin_counts = ni)

unb_row <- data.frame(as.character(round(lambda0, 4)),
                      unb_res[1],unb_res[2],unb_res[3],unb_res[4])

######## Creating results table ################################

colnames(res_sig_search) <- colnames(sf_row) <- colnames(unb_row) <- 
  c("$\\lambda$ / Test", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif',
    '$\\hat{S}$')

res_sig_search <- rbind(res_sig_search, unb_row, sf_row)

delta_tab_res <- sapply(c(lambda_seq, lambda0), function(t){
  unbinned_delta_estimation(lambda = t, fs = fs, qb = qb,
                            mu1 = mean1_in_gb, mu2 = mean2_in_gb,
                            sd = sd_in_gb, search_region = c(l,u),
                            bkg_data = bkg_only_data)
}) |> cbind(c(delta_hat_sf, delta_se_sf))|> round(4)

delta_col <- apply(delta_tab_res, 2, function(x) {paste0(x[1], " +/- " , x[2])})
res_sig_search <- cbind(res_sig_search, delta_col)
res_sig_search <- res_sig_search[,
                                 c(1,ncol(res_sig_search),2:(ncol(res_sig_search)-1))] # rearranging the columns
colnames(res_sig_search)[2] <- "$\\hat{\\delta}$"

caption <- paste0("Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2, '<br> Bkg: BW X Power')
kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center')

######### Plotting Densities ###################################
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
        col = alpha(palette()[j], alpha = 0.6),
        lty = my_lty[j])
}

curve(N*gb_sf(x)*(u-l)/(k+1),
      l, u, add = TRUE, lwd = 2.2,
      col = alpha('blue', alpha = 0.6),
      lty = 1)

curve(N*((dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean1_in_gb, sd = sd_in_gb) + 
            dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mean2_in_gb, sd = sd_in_gb))*lambda0 +
           (1-2*lambda0)*qb(x))*(u-l)/(k+1),
      l, u, add = TRUE, lwd = 2.2,
      col = alpha('green', alpha = 0.6),
      lty = 1)

legend('bottomleft', col = c('blue', mycols, 'green'),
       lty = c(1:6,1), bty = 'n', lwd =2.2,
       legend=c('safeguard',
                TeX(sprintf(r'($g_b(\lambda = %f)$)', c(lambda_seq, 
                                                        round(lambda0, 4))))),
       cex = 1,
       y.intersp = 1)

