rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

mean_sig <- 125; sd_sig <- 2
eps <- 1e-3
u <- 160; l <- 110
n_breaks = 50
obs <- read.csv('../Data/toydata_sig10_bkg500.csv', header = TRUE)$x

hist(obs, probability = TRUE, breaks = n_breaks)

n <- length(obs)
obs_normalised <- scale(obs)
kde_ <- kdensity(obs_normalised, bw = 0.1)
kde_unnormed<- function(t) kde_((t-mean(obs))/sd(obs))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
curve(kde, col = 'blue', add = TRUE, lwd = 2)

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



bkg_loc <- 91.2
bkg_scale <- 2.49/2

qb <- function(x)
{
  dtrunc(x, a = l, b= u, spec = 'cauchy', 
         location = bkg_loc, 
         scale = bkg_scale)
}

hist(obs, probability = TRUE, breaks = n_breaks)
curve(kde, lwd = 2, col = 'blue', add = TRUE)
curve(qb, lwd = 2, col = 'brown', add = TRUE)

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- qb(x)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

# Area under gb
integrate(gb_test,l,u)

search_signal <- function(lambda)
{
  norm_S <- integrate(function(t) {((fs(t)/gb_test(t, fs_prop = lambda) - 1)^2)*gb_test(t, fs_prop = lambda)}, l, u)$value |> sqrt()
  
  S1_vec <- sapply(obs, function(x) {
    f_sig <- fs(x)
    g_b <- gb_test(x, fs_prop = lambda)
    return((f_sig/g_b -1)/norm_S)
  })
  theta <- mean(S1_vec)
  se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)
  
  # testing \eta = 0:
  eta_hat <- theta/norm_S
  
  # testing for signal:
  (t_stat_theta <- theta/se_theta)
  p_val <- pnorm(t_stat_theta, lower.tail = FALSE)
  
  return(c(eta_hat, p_val))
  
}

lambda_max <- 0.005
(lambda_seq <- seq(0, lambda_max, length.out = 4))
res_sig_search <- sapply(lambda_seq, search_signal)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
sig_lev <- gtools::stars.pval(res_sig_search[,3])
res_sig_search <- cbind(res_sig_search, sig_lev)
colnames(res_sig_search) <- c('lambda', 'eta_hat', 'p-value', '')
res_sig_search


mycols <- c('red', 'green', 'orange', 'purple')
palette(mycols)
my_lty = c(1,2,4,5)

curve(kde, M_lower - 1, M_upper + 1, xlab = 'y',
      ylab = 'Density', lwd = 2.2, ylim = c(kde(M_upper+1), kde(M_lower-1)))

for(j in 1:length(lambda_seq))
{
  curve(gb_test(x, fs_prop = lambda_seq[j]),
        M_lower - 1, M_upper + 1, add = TRUE, lwd = 2.2,
        col = j,
        lty = my_lty[j])
}

abline(v = c(M_lower, M_upper), col = 'black', lwd = 2)

legend('bottomleft', col = 1:length(lambda_seq),
       lty = my_lty, bty = 'n', lwd =2.2, 
       legend=TeX(sprintf(r'($\lambda = %f (%f)$)', lambda_seq, 
                          round(as.numeric(res_sig_search[,3]), 4))),
       cex = 1.5)

#####################################################################
######################## ADDING SAFEGUARD RESULTS ###################
#####################################################################
set.seed(1008)

bkg_sample <- rtrunc(n, a = l, b = u,
                     spec = 'cauchy',
                     location = bkg_loc,
                     scale = bkg_scale)

bkg_mod <- function(pars, dat, mean = mean_sig)
{
  lam <- pars
  fi <- lam*dtrunc(dat, a = l, b = u, spec = 'norm',
                   mean = mean, sd = sd_sig) + 
    (1-lam)*dtrunc(dat,a=l,b=u,spec="cauchy",
                   location = bkg_loc,
                   scale = bkg_scale)
  return(-sum(log(fi)))
}

bounds <- matrix(c(
  0, 1
), ncol = 2, byrow = TRUE)

colnames(bounds) <- c("lower", "upper")
K <- nrow(bounds)
ui <- rbind( diag(K), -diag(K) )
ci <- c( bounds[,1], - bounds[,2] )
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

mean_seq <- c(seq(M_lower, mean_sig, length.out = 25),
              seq(mean_sig, M_upper, length.out = 25)[-1])
res_list <- sapply(mean_seq, function(x) nlminb(start = 0.1,
                                                objective = bkg_mod,
                                                upper = 1, lower = 0,
                                                mean = x, dat = bkg_sample)$par)

(lam_hat <- max(res_list))

gb_sf <- function(x){
  lam_hat*dtrunc(x, a = l, b = u, spec = 'norm',
                 mean = mean_sig, sd = sd_sig) + 
    (1-lam_hat)*dtrunc(x, a = l, b = u, spec = 'cauchy',
                       location = bkg_loc, scale = bkg_scale)
}

mix_mod_safeguard <- function(pars, dat)
{
  eta <- pars
  mix_dens <- log(eta*dtrunc(dat, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig) + 
                    (1-eta)*(lam_hat*dtrunc(dat, a = l, b = u, spec = 'norm', mean = mean_sig, sd = sd_sig) + 
                               (1-lam_hat)*dtrunc(dat,a=l,b=u,spec="cauchy",
                                                  location = bkg_loc,
                                                  scale = bkg_scale)))
  log_lik <- sum(mix_dens)
  return(-log_lik)
}

res_safeguard <- nlminb(start = 0.1,
                        objective = mix_mod_safeguard,
                        lower = 0,
                        upper = 1,
                        dat = obs)
(eta_safeguard <- res_safeguard$par)
ll0 <- -mix_mod_safeguard(pars = 0, dat = obs)
ll1 <- -res_safeguard$objective
test_stat_safeguard <- -2*(ll0-ll1)
(p_val_safeguard <- 0.5*pchisq(test_stat_safeguard, df = 1, lower.tail = FALSE))


hs <- hist(obs, probability = TRUE, breaks = n_breaks)
curve(kde, lwd = 2, col = 'blue', add = TRUE)
curve(qb, lwd = 2, col = 'red', add = TRUE)
curve(gb_sf, lwd = 2, col = 'cyan', add = TRUE)
abline(v = c(M_lower, M_upper), col = 'black', lwd = 2)

################################################################################
###################### Making the Binned plot ##################################
################################################################################

picture_l <- M_lower - 3; picture_u <- M_upper + 3
plot(hs$mids, hs$counts, pch = 16,
     ylab="Events",xlab="Energy (Gev)",
     xlim = c(picture_l, picture_u))
arrows(x0 = hs$mids, y0 = hs$counts + sqrt(hs$counts),
       x1 = hs$mids, y1 = hs$counts - sqrt(hs$counts),
       angle = 90, length = 0.01, code = 3, col = 'red')

curve(n*kde(x)*(u-l)/n_breaks, lwd = 2.2, col = 'black', add = TRUE, lty = 1)
curve(n*gb_sf(x)*(u-l)/n_breaks, lwd = 2.2, col = 'cyan', add = TRUE, lty = 2)
abline(v = c(M_lower, M_upper), col = 'black', lwd = 2.2)

mycols <- c('red', 'green', 'orange', 'purple')
palette(mycols)
my_lty = c(3,4,5,6)

for(j in 1:length(lambda_seq))
{
  curve(n*gb_test(x, fs_prop = lambda_seq[j])*(u-l)/n_breaks,
        l, u, add = TRUE, lwd = 2.2,
        col = j,
        lty = my_lty[j])
}

legend('topright', col = c('black', 'cyan', mycols),
       lty = c(1, 2, my_lty), bty = 'n', lwd =2.2,
       legend=c('KDE',
                TeX(sprintf(r'(safeguard $(\hat{\eta}: %f, p:%f)$)', eta_safeguard, p_val_safeguard)),
                TeX(sprintf(r'($\lambda = %f (\hat{\eta}: %f , p : %f)$)', lambda_seq,
                            round(as.numeric(res_sig_search[,2]), 4),
                            round(as.numeric(res_sig_search[,3]), 4)))),
       cex = 1,
       y.intersp = 1)
