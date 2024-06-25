rm(list = ls())
library(ggplot2)
source('functions truncnorm without eta.R')
# data generation
set.seed(12345)
n <- 2e3
back <- rtruncnorm(n, a = l, b= u, mean = mean_back,
                   sd = sd_back)
sig <- rtruncnorm(n, a = l, b= u, mean = mean_sig,
                  sd = sd_sig)
u_mask <- runif(n)
obs <- ifelse(u_mask<1-eta_true,back,sig)

xs <- seq(l,u,0.01)
y_sig <- dtruncnorm(xs, mean = mean_sig, sd = sd_sig, 
                    a = l, b = u)
y_back <- dtruncnorm(xs,mean = mean_back, sd = sd_back, 
                     a = l, b = u)
y_real <- eta_true*y_sig + (1 - eta_true)*y_back


hs <- hist(obs, breaks = 20, probability = TRUE)

ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                  aes(y = after_stat(density)),
                                                fill = 'steelblue', col = 'black',
                                                breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = 1, color = 'Proposal Background'),
                                   lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black'))

bias_vector <- c()
# f = g_b[1 + S1 + T1]-----------------------------------------------
t <- Sys.time()
res1 <- constrOptim(theta = rep(0.01,2),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0), c(-1,0)),
                    ci = c(0,-norm_S))

(theta_1 <- res1$par[1])
(tau_1 <- res1$par[-1])

(eta_1_biased <- theta_1/norm_S)

(bias_1_observed <- eta_1_biased - eta_true)

# (delta_1 <- (theta_1 - eta_true*norm_S)/(1-eta_true))
theta_MME <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                          sd = sd_sig)) |> mean()
(delta_MME <- (theta_MME - eta_true*norm_S)/(1-eta_true))
(delta_true <- sapply(back, function(t) S1(t, mean = mean_sig, sd = sd_sig))|> mean())
(bias_theoretical <- (1-eta_true)*delta_true/norm_S)
bias_vector <- c(bias_vector, bias_1_observed)
(time_taken_1 <- Sys.time() - t)

# f = g_b[1 + S1 + T1 + T2]---------------------------------------------------
t <- Sys.time()
res2 <- constrOptim(theta = rep(0.01,3),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0), c(-1,0,0)),
                    ci = c(0,-norm_S))

(theta_2 <- res2$par[1])
(tau_2 <- res2$par[-1])

(eta_2_biased <- theta_2/norm_S)

(bias_2_observed <- eta_2_biased - eta_true)
bias_vector <- c(bias_vector, bias_2_observed)
bias_theoretical

(time_taken_2 <- Sys.time() - t)

# f = g_b[1 + S1 + T1 + T2 + T3]----------------------------------------------
t <- Sys.time()
res3 <- constrOptim(theta = rep(0.01,4),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0), c(-1,0,0,0)),
                    ci = c(0,-norm_S))

(theta_3 <- res3$par[1])
(tau_3 <- res3$par[-1])

(eta_3_biased <- theta_3/norm_S)

(bias_3_observed <- eta_3_biased - eta_true)
bias_vector <- c(bias_vector, bias_3_observed)
bias_theoretical

(time_taken_3 <- Sys.time() - t)

# f = g_b[1 + S1 + T1 + T2 + T3 + T4]----------------------------------------------
# t <- Sys.time()
# res4 <- constrOptim(theta = rep(0.1,5),f = neg_loglikelihood,
#                     data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
#                     grad = score, 
#                     ui = rbind(c(1,0,0,0,0), c(-1,0,0,0,0)),
#                     ci = c(0,-norm_S))
# 
# # Error in optim(theta.old, fun, gradient, control = control, method = method,  : 
# # initial value in 'vmmin' is not finite
# 
# (theta_4 <- res4$par[1])
# (tau_4 <- res4$par[-1])
# 
# (eta_4_biased <- theta_4/norm_S)
# 
# (bias_4_observed <- eta_4_biased - eta_true)
# bias_theoretical
# 
# (time_taken_4 <- Sys.time() - t)

#-----------------------------------------------------------------------------

plot(bias_vector, type = 'o', ylim = c(-0.006,-0.0048))
abline(h = bias_theoretical, col = 'red', lty = 2)

save.image(file = "fitting truncnorm without eta.RData")
