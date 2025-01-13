rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

mean_sig <- 125; sd_sig <- 2
eps <- 1e-3
u <- 160; l <- 110

obs <- read.csv('Data/toy3.csv', header = TRUE)$x

hist(obs, probability = TRUE, breaks = 50)

n <- length(obs)
obs_normalised <- scale(obs)
kde_ <- kdensity(obs_normalised)
kde_unnormed<- function(t) kde_((t-mean(obs))/sd(obs))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot

curve(kde, lwd = 2, col = 'blue', add = TRUE)

# qb_ <- function(t) dtrunc(t, spec = 'exp', rate = 0.05, a = 0, b = u-l)
# 
# qb <- function(t) qb_(t-l)
# 
# curve(qb, lwd = 2, col = 'brown', add = TRUE)

qb <- function(t) dtrunc(t, spec = 'pareto', shape = 7, scale = l,
                         a = l, b = u)

hist(obs, probability = TRUE, breaks = 50)
curve(kde, lwd = 2, col = 'blue', add = TRUE)
curve(qb, lwd = 2, col = 'brown', add = TRUE)

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

M_lower <- mean_sig - r
M_upper <- mean_sig + r

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm')
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm')
  qb_val <- qb(x)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

# Area under gb
integrate(gb_test,l,u)

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 0.1) +
  stat_function(fun = kde, col = 'orange') +
  stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = function(t) gb_test(t), col = 'red') + 
  geom_vline(xintercept = c(M_lower, M_upper))


mycols <- c('red', 'blue', 'green', 'brown', 'orange')
palette(mycols)
my_lty = c(1,2,4,5,3)
fs_prop_seq <- c(0, seq(0.01, 0.04, 0.01))
curve(kde, M_lower - 5, M_upper + 5, xlab = 'y',
      ylab = 'Density', lwd = 2.2)
for(j in 1:length(fs_prop_seq))
{temp_fun <- function(t) gb_test(t, fs_prop = fs_prop_seq[j])
curve(temp_fun, M_lower - 5, M_upper + 5, add = TRUE, lwd = 2.2,
      col = j,
      lty = my_lty[j])
}

legend(x = 130, y = 0.04, col = 1:length(fs_prop_seq),
       lty = my_lty, bty = 'n', lwd =2.2, 
       legend=TeX(sprintf(r'($\lambda = %f$)', fs_prop_seq)),
       cex = 1.5)

gb <- function(t) gb_test(t, fs_prop = 0.01)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 0.1) +
  stat_function(fun = kde, col = 'black', lwd = 1.5) + 
  # stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = gb, col = 'orange', lwd = 1.5) + 
  geom_vline(xintercept = c(M_lower, M_upper), alpha = 0.3,
             lwd = 1.5) + 
  xlab('x') + ylab('Density') +
  annotate('text', x = c(M_lower-0.1, M_upper+0.1),
           y = c(0,0), 
           label = c(TeX('$\\mu_s-d_\\epsilon$'), TeX('$\\mu_s+d_\\epsilon$')),
           size = 10) + 
  theme_bw() + My_Theme

# Function to calculate norm under Gb
calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

# Convarting fs into S1:
norm_S <- calc_norm_gb(function(t) fs(t)/gb(t) - 1)

S1 <- function(x, ...)
{
  f_sig <- fs(x, ...)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

S1_vec <- sapply(obs, S1)
theta <- mean(S1_vec)

se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)

# testing \eta = 0:
theta/norm_S

# testing for signal:
(t_stat_theta <- theta/se_theta)
pnorm(t_stat_theta, lower.tail = FALSE)
t_stat_theta>qnorm(0.95)
