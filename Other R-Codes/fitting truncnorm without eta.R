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
# f = g_b[1 + S1 + T1]---------------------------------------------------
t <- Sys.time()
res1 <- constrOptim(theta = rep(0.01,2),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0), c(-1,0)),
                    ci = c(0,-norm_S))

(time_taken_1 <- Sys.time() - t)
(theta_1 <- res1$par[1])
(tau_1 <- res1$par[-1])
(eta_1 <- theta_1/norm_S)

fit_1 <- sapply(xs, function(t) mod(t, theta = theta_1, tau = tau_1,
                                    mean = mean_sig, sd = sd_sig))
mix_1 <- as.numeric(fit_1[1,])
back_from_mix_1 <- as.numeric(fit_1[2,])

plt_1 <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                               aes(y = after_stat(density)),
                                                             fill = 'steelblue', col = 'black',
                                                             breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = 1, color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = mix_1, color = 'Fitted Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = back_from_mix_1, color = 'Fitted Background from Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black',
                                'Fitted Mixture' = 'cyan',
                                'Fitted Background from Mixture' = 'skyblue'))+ 
  ggtitle(bquote(T[1]))

plt_1

# f = g_b[1 + S1 + T1 + T2]---------------------------------------------------
t <- Sys.time()
res2 <- constrOptim(theta = rep(0.01,3),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0), c(-1,0,0)),
                    ci = c(0,-norm_S))

(time_taken_2 <- Sys.time() - t)
(theta_2 <- res2$par[1])
(tau_2 <- res2$par[-1])
(eta_2 <- theta_2/norm_S)

fit_2 <- sapply(xs, function(t) mod(t, theta = theta_2, tau = tau_2,
                                    mean = mean_sig, sd = sd_sig))
mix_2 <- as.numeric(fit_2[1,])
back_from_mix_2 <- as.numeric(fit_2[2,])

plt_2 <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                               aes(y = after_stat(density)),
                                                             fill = 'steelblue', col = 'black',
                                                             breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = 1, color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = mix_2, color = 'Fitted Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = back_from_mix_2, color = 'Fitted Background from Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black',
                                'Fitted Mixture' = 'cyan',
                                'Fitted Background from Mixture' = 'skyblue'))+ 
  ggtitle(bquote(T[1] + T[2]))

plt_2

# f = g_b[1 + S1 + T1 + T2 + T3]----------------------------------------------
t <- Sys.time()
res3 <- constrOptim(theta = rep(0.01,4),f = neg_loglikelihood,
                    data = obs, mean_sig = mean_sig, sd_sig = sd_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0), c(-1,0,0,0)),
                    ci = c(0,-norm_S))

(time_taken_3 <- Sys.time() - t)
(theta_3 <- res3$par[1])
(tau_3 <- res3$par[-1])
(eta_3 <- theta_3/norm_S)

fit_3 <- sapply(xs, function(t) mod(t, theta = theta_3, tau = tau_3,
                                    mean = mean_sig, sd = sd_sig))
mix_3 <- as.numeric(fit_3[1,])
back_from_mix_3 <- as.numeric(fit_3[2,])

plt_3 <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                               aes(y = after_stat(density)),
                                                             fill = 'steelblue', col = 'black',
                                                             breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = 1, color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = mix_3, color = 'Fitted Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = back_from_mix_3, color = 'Fitted Background from Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black',
                                'Fitted Mixture' = 'cyan',
                                'Fitted Background from Mixture' = 'skyblue'))+ 
  ggtitle(bquote(T[1] + T[2] + T[3]))

plt_3
#-----------------------------------------------------------------------------

final_plt <- ggpubr::ggarrange(plt_1, plt_2, plt_3,
                                   ncol = 3, nrow = 1,
                                   common.legend = TRUE, legend = 'bottom')
final_plt

save.image(file = "fitting truncnorm without eta.RData")
