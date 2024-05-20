library(truncdist)
library(VGAM)
library(ggplot2)

source('bases on [l,u] with tuned gb.R')

# BACKGROUND AND MIXTURE PARAMETERS:
shape_back <- 1.4; scale_back <- 1
eta_true <- 0

fb <- function(x)
{
  dtrunc(exp(x), shape = shape_back,
         scale = scale_back, spec = 'pareto', 
         a = real_l, b = real_u)*exp(x)
}

# checking delta:
(true_delta <- integrate(function(t) S1(t)*fb(t), l, u)$value)


N_rep <- 1e3
n <- 2338

theta_vec <- c()
t_vec <- c()

data_list <- list()

set.seed(1234)

time_now <- Sys.time()
for(i in 1:N_rep)
{
  samp_sim <- rtrunc(n, mean = mean_sig, sd = sd_sig,
                     spec = 'norm', a = real_l, b = real_u)
  bkg_sim <- rtrunc(n, shape = shape_back,
                    scale = scale_back, spec = 'pareto',
                    a = real_l, b = real_u)
  u_mask <- runif(n)
  
  obs_sim <- ifelse(u_mask<eta_true, samp_sim, bkg_sim)
  obs_sim <- log(obs_sim)
  data_list[[i]] <- obs_sim
  
  S1_vec <- sapply(obs_sim, function(t) S1(t))
  theta <- mean(S1_vec)
  se <- sqrt((mean(S1_vec^2) - theta^2)/n)
  
  theta_vec <- c(theta_vec, theta)
  t_vec <- c(t_vec, theta/se)
  
  message(sprintf('Iteration: %d/%d',i,N_rep))
}

Sys.time() - time_now

mean(t_vec > qnorm(0.95))

mean(theta_vec<0)

ggplot(mapping = aes(y = theta_vec/norm_S, x = 1:N_rep)) + 
  geom_point() + geom_line()

hist(t_vec, probability = TRUE)

hist(theta_vec, probability = TRUE)
abline(v = true_delta, col = 'red')
