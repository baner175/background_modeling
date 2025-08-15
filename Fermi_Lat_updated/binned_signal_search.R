rm(list = ls())

library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)

real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
eps <- 1e-3

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

# SIGNAL CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd))
}

dat <- read.table('Data_ex1.txt', header = TRUE)
x <- dat$x
y <- log(x)
N <- length(y)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
ni <- sapply(1:k, function(i){
  sum(y>bin_ends[i] & y<=bin_ends[i+1])
})

qb_y_model <- function(beta){
  qb_i <- sapply(1:k, function(i){
    integrate(function(x){
      dtrunc(x, spec = 'exp', rate = beta, a = l, b = u)
    }, bin_ends[i], bin_ends[i+1])$value
  })# pareto model on x is the same as exponential model on y = log(x)
  return(-sum(ni*log(qb_i)))
}

beta_hat <- nlminb(start = 0.01,
                   objective = qb_y_model,
                   upper = Inf, lower = -Inf)$par

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(log(mean_sig)-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(log(mean_sig)+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(log(mean_sig) - l,u - log(mean_sig)))

r <- sol$root

M_lower <- log(mean_sig) - r
M_upper <- log(mean_sig) + r

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps


plot(x = xi, y = ni, pch = 16, col = 'grey',
     xlab = 'log(x)', ylab = 'counts')
curve(N*dtrunc(x, spec = 'exp', rate = beta_hat, a = l, b = u)*(u-l)/(k+1),
      lwd = 2, col = 'blue', add = TRUE)
curve(N*fs(x)*(u-l)/(k+1), add = TRUE, col = 'green')
abline(v = c(M_lower, M_upper), col = 'black', lty = 2, lwd = 2)

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (exp(M_lower) + mean_sig)/2
mean2_in_gb <- (exp(M_upper) + mean_sig)/2
sd_in_gb <- 2*sd_sig

signal_search <- function(lambda){
  norm_S <- integrate(function(y){
    fs <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                 mean = mean_sig, sd = sd_sig)*exp(y)
    qb <- dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
    fs_val1 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean1_in_gb, sd = sd_in_gb)*exp(y)
    fs_val2 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean2_in_gb, sd = sd_in_gb)*exp(y)
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    return(((fs/gb - 1)^2)*gb)
  }, l, u)$value |> sqrt()
  
  S2_vec <- sapply(xi, function(y){
    fs <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                 mean = mean_sig, sd = sd_sig)*exp(y)
    qb <- dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
    fs_val1 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean1_in_gb, sd = sd_in_gb)*exp(y)
    fs_val2 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean2_in_gb, sd = sd_in_gb)*exp(y)
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    S_val <- (fs/gb - 1)
    return(S_val/(norm_S^2))
  })
  
  theta0_hat <- sum(S2_vec*ni)/N
  
  d_log_qb_xi <- (1/beta_hat) - xi - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  
  d2_log_qb_xi <- -(1/(beta_hat^2)) - (((l^2)*exp(-beta_hat*l) - (u^2)*exp(-beta_hat*u))/(exp(-beta_hat*l) - exp(-beta_hat*u)) - 
                                        ((u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)))^2)
  
  d_normS_sq <- -(1-2*lambda)*integrate(function(y){
    fs <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                 mean = mean_sig, sd = sd_sig)*exp(y)
    qb <- dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
    fs_val1 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean1_in_gb, sd = sd_in_gb)*exp(y)
    fs_val2 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean2_in_gb, sd = sd_in_gb)*exp(y)
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- (1/beta_hat) - y - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
    return(((fs/gb)^2)*qb*d_log_qb)
  }, l, u)$value
  
  d_S2_xi <- sapply(xi, function(y){
    fs <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                 mean = mean_sig, sd = sd_sig)*exp(y)
    qb <- dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
    fs_val1 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean1_in_gb, sd = sd_in_gb)*exp(y)
    fs_val2 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                      mean = mean2_in_gb, sd = sd_in_gb)*exp(y)
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- (1/beta_hat) - y - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
    
    return(-((norm_S^2)*(fs/(gb^2))*(1-2*lambda)*qb*d_log_qb + (fs/gb-1)*d_normS_sq)/(norm_S^4))
  })
  
  J_hat <- (-1/k)*sum(ni*d2_log_qb_xi)
  d_theta0 <- sum(ni*d_S2_xi)/N
  var_S2_F_hat <- sum((S2_vec^2)*ni)/N - theta0_hat^2
  c_hat <- N/k
  
  
  sig_theta0_hat <- sqrt(
    var_S2_F_hat + 
      ((c_hat^2)/(J_hat^2))*(d_theta0^2)*(sum(ni*(d_log_qb_xi)^2)/N) + 
      (2*c_hat/J_hat)*d_theta0*
      (sum(S2_vec*d_log_qb_xi*ni)/N)
  )
  
  theta0_stat <- sqrt(N)*(theta0_hat-0)/sig_theta0_hat
  p_val <- pnorm(theta0_stat, lower.tail = FALSE)
  # S_hat <- n*theta0_hat
  
  return(c(theta0_hat,
           p_val,
           qnorm(p_val, lower.tail = FALSE) |> round(1)
  ))
}

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 6)

res_sig_search <- sapply(lambda_seq, signal_search)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c("$\\lambda$", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif')
caption <- paste0("Signal Search Results (Binned Data)")
kable(res_sig_search, "html", digits = 5,
      booktabs = TRUE, escape = FALSE,
      caption = caption) |>
  kable_styling(full_width = FALSE, position = 'center')

######### Plotting Densities ###################################################

qb <- function(y) dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
gb <- function(y, lambda){
  fs <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
               mean = mean_sig, sd = sd_sig)*exp(y)
  qb <- dtrunc(y, spec = 'exp', rate = beta_hat, a = l, b = u)
  fs_val1 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                    mean = mean1_in_gb, sd = sd_in_gb)*exp(y)
  fs_val2 <- dtrunc(exp(y), spec = 'norm', a = real_l, b = real_u,
                    mean = mean2_in_gb, sd = sd_in_gb)*exp(y)
  gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
  return(gb)
}

plot(y = ni, x = xi,
     pch = 16, col = 'grey', xlab = 'log(x)',
     main = TeX('Sensitivity Analysis for $g_b(x; \\hat{beta}_{bkg}, \\lambda)$'))

picture_l <- M_lower - 3; picture_u <- M_upper + 3

mycols <- c('black','brown', 'skyblue' ,'red', 'orange', 'purple')
palette(mycols)
my_lty = 1:6
for(j in 1:length(lambda_seq))
{
  curve(N*gb(x, lambda = lambda_seq[j])*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 2.2,
        col = ggplot2::alpha(palette()[j], alpha = 0.6),
        lty = my_lty[j])
}

legend('topright', col = mycols,
       lty = 1:5, bty = 'n', lwd =2.2,
       legend=TeX(sprintf(r'($g_b(\lambda = %f)$)', lambda_seq)),
       cex = 1,
       y.intersp = 1)


