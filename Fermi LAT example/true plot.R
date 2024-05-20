rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)

real_l <- 1; real_u <- 35

l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
phi_bkg <- 1.4
eta_true <- 0.02


dat <- read.table('Data_ex1.txt', header = TRUE)
obs <- log(dat$x)
n <- length(obs)

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

# BACKGROUND DENSITY
fb <- function(x)
{
  return(dtrunc(exp(x), spec = 'pareto', 
                shape = phi_bkg, scale = real_l,
                a = real_l, b = real_u)*exp(x))
}

# MISTURE DENSITY
f_mix <- function(x)
{
  return((1-eta_true)*fb(x)+eta_true*fs(x))
}

dat <- read.table('Data_ex1.txt', header = TRUE)
obs <- log(dat$x)
n <- length(obs)
xs <- seq(l, u, length.out = 1e3)
y_bkg <- sapply(xs, fb)
y_mix <- sapply(xs, f_mix)

hs <- hist(obs, probability = TRUE, breaks = 50)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

basic_plt <- ggplot(mapping = aes(x = obs)) + 
  geom_histogram(aes(y = after_stat(density)),
                 fill = 'grey', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(aes(x = xs, y = y_bkg), col = 'red', linetype = 2,
            lwd = 1.2) +
  geom_line(aes(x = xs, y = y_mix), col = 'blue', linetype = 1,
            lwd = 1.2) + xlab('log(y)')+ylab('Density') +
  theme_bw() + 
  My_Theme
basic_plt



