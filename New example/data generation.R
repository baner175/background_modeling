rm(list = ls())
library(truncdist)
library(ggplot2)
set.seed(1234567)

l <- 1; u <- 5
mean_sig <- 2.5; sd_sig <- 0.1
mean_back <- 0.5; sd_back <- 2.5
eta_true <- 0.03
n <- 5e3
eps <- 1e-3

# Use larger sample
# Do a 50-50 split
# do a stratified split

sig <- rtrunc(n, spec = 'norm', a = l, b = u,
              mean = mean_sig, sd = sd_sig)
back <- rtrunc(n, spec = 'norm', a = l, b = u,
               mean = mean_back, sd = sd_back)
u_mask <- runif(n)
obs <- ifelse(u_mask<eta_true, sig, back)

hs <- hist(obs, breaks = 20, probability = TRUE)

fs <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                            mean = mean_sig, sd = sd_sig)
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

f_back <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                             mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*f_back(t)

xs <- seq(l, u, 0.01)
ggplot(mapping = aes(x = obs)) + 
  geom_histogram(mapping = aes(y = after_stat(density)), breaks = hs$breaks,
                 fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = xs, y = f_back(xs)), col = 'brown',
            lwd = 1.2, alpha = 0.5) + 
  geom_line(mapping = aes(x = xs, y = f_mix(xs)), col = 'black',
            lwd = 1.2, alpha = 0.5)

find_d <- function(d)
{
  pl <- Fs(mean_sig-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(mean_sig+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

M_lower <- mean_sig - r
M_upper <- mean_sig + r

sig_mask <- obs<=M_upper & obs>=M_lower
data_sig <- obs[sig_mask]
data_bkg <- obs[!sig_mask]

sig_train_index <- sample(1:length(data_sig), replace = FALSE, size = 0.5*length(data_sig))
bkg_train_index <- sample(1:length(data_bkg), replace = FALSE, size = 0.5*length(data_bkg))

train_sig <- data_sig[sig_train_index]
train_bkg <- data_bkg[bkg_train_index]
train <- c(train_sig, train_bkg)

test_sig <- data_sig[-sig_train_index]
test_bkg <- data_bkg[-bkg_train_index]
test <- c(test_sig, test_bkg)


hist(obs, probability = TRUE, breaks = 50)
curve(f_mix, add = TRUE)
hist(test, probability = TRUE, breaks = 50)
curve(f_mix, add = TRUE)
hist(train, probability = TRUE, breaks = 50)
curve(f_mix, add = TRUE)

# write.csv(data.frame(x = obs), 'full_data.csv', row.names = FALSE)
# write.csv(data.frame(x = train), 'training_data.csv', row.names = FALSE)
# write.csv(data.frame(x = test), 'test_data.csv', row.names = FALSE)
  



