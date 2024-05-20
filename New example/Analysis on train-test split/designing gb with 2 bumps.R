rm(list = ls())

library(truncdist)
library(ggplot2)
obs <- read.csv('training_data.csv', header = TRUE)$x
n <- length(obs)

l <- 1; u <- 5
mean_sig <- 2.5; sd_sig <- 0.1
rate_gb <- 0.35
eps <- 1e-3

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

find_d <- function(d)
{
  pl <- Fs(mean_sig-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(mean_sig+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

(M_lower <- mean_sig - r)
(M_upper <- mean_sig + r)

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps

hs <- hist(obs, probability = TRUE, breaks = 50)
ggplot(mapping = aes(x = obs)) + 
  geom_histogram(aes(y = after_stat(..density..)),
                 breaks = hs$breaks,
                 fill = 'steelblue', col = 'black')

kde <- kdensity::kdensity(obs)

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 0.6) +
  stat_function(fun = kde, col = 'orange') + 
  stat_function(fun = fs, col = 'blue') + 
  geom_vline(xintercept = c(M_lower, M_upper))

mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 1.9*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val_1 <- fs(x, mean = mean1_in_gb, sd = sd_in_gb)
  fs_val_2 <- fs(x, mean = mean2_in_gb, sd = sd_in_gb)
  gb_val <- dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u)
  return((fs_prop)*fs_val_1 + (fs_prop)*fs_val_2 + (1-2*fs_prop)*gb_val)
}

integrate(gb_test, l, u)

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 0.5) +
  stat_function(fun = kde, col = 'orange') + 
  stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = function(t) gb_test(t), col = 'red') + 
  geom_vline(xintercept = c(M_lower, M_upper))

curve(kde, l, u)
abline(v = c(M_lower, M_upper))
fs_prop_seq <- seq(0, 0.03, 0.005)
for(j in 1:length(fs_prop_seq))
{temp_fun <- function(t) gb_test(t, fs_prop = fs_prop_seq[j])
curve(temp_fun, l, u, add = TRUE, col = j, lty = 2)
}

legend('topright', col = 1:length(fs_prop_seq),
       lty = 2, bty = 'n', legend = as.character(fs_prop_seq))

# Probably best to use 0.02

gb <- function(t) gb_test(t, fs_prop = 0.02)

ggplot(data=data.frame(x = c(l,u)), mapping = aes(x = x)) + 
  ylim(0, 0.5) +
  stat_function(fun = kde, col = 'orange') + 
  stat_function(fun = fs, col = 'blue') + 
  stat_function(fun = gb, col = 'red') + 
  geom_vline(xintercept = c(M_lower, M_upper))


