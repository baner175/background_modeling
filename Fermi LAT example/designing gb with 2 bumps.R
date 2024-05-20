rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

real_l <- 1; real_u <- 35

l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
rate_gb <- 1.2
eps <- 1e-3


dat <- read.table('Data_ex1.txt', header = TRUE)
obs <- log(dat$x)
n <- length(obs)
kde <- kdensity(obs)

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

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (exp(M_lower) + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (exp(M_upper) + mean_sig)/2;
gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- fs(x, mean = mean1_in_gb, sd = sd_in_gb)
  fs_val2 <- fs(x, mean = mean2_in_gb, sd = sd_in_gb)
  qb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

# Area under gb
integrate(gb_test,l,u)

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 1.8) +
  stat_function(fun = kde, col = 'orange') + 
  stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = function(t) gb_test(t), col = 'red') + 
  geom_vline(xintercept = c(M_lower, M_upper))


mycols <- c('red', 'blue', 'green', 'brown', 'orange')
palette(mycols)
my_lty = c(1,2,4,5,3)
fs_prop_seq <- seq(0, 0.04, 0.01)
curve(kde, l, u, xlab = 'log(y)',
      ylab = 'Density', lwd = 2.2)
for(j in 1:length(fs_prop_seq))
{temp_fun <- function(t) gb_test(t, fs_prop = fs_prop_seq[j])
  curve(temp_fun, l, u, add = TRUE, col = j, lty = my_lty[j], lwd = 2.2)
}

legend(x = 2, y = 0.9, col = 1:length(fs_prop_seq),
       lty = my_lty, bty = 'n', lwd =2.2, 
       legend=TeX(sprintf(r'($\lambda = %f$)', fs_prop_seq)),
       cex = 1.5)


# Probably best to use 0.02

gb <- function(t) gb_test(t, fs_prop = 0.02)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 1.8) +
  stat_function(fun = kde, col = 'black', lwd = 1.5) + 
  # stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = gb, col = 'orange', lwd = 1.5) + 
  geom_vline(xintercept = c(M_lower, M_upper), alpha = 0.5,
             lwd = 1.5) + 
  xlab('x') + ylab('Density') +
  annotate('text', x = c(M_lower, M_upper),
           y = c(0,0), 
           label = c(TeX('$\\mu_s-d_\\epsilon$'), TeX('$\\mu_s+d_\\epsilon$')),
           size = 10) + 
  theme_bw() + My_Theme

