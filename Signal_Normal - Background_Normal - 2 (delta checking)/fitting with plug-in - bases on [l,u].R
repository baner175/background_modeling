# source('bases on [l,u].R')
source('bases on [l,u] with tuned gb.R')

library(ggplot2)

# checking if delta < 0:
integrate(function(t) S1(t, mean = mean_sig,
                         sd = sd_sig)*dtrunc(t, mean = mean_back,
                                             sd = sd_back, spec = 'norm',
                                             a = l, b = u),
          lower = l, upper = u)$value

xs <- seq(l,u,0.01)
y_sig <- dtrunc(xs, spec = 'norm',
                mean = mean_sig, sd = sd_sig,
                a = l, b = u)
y_back <- dtrunc(xs, spec = 'norm', 
                 mean = mean_back, sd = sd_back,
                 a = l, b = u)
y_real <- eta_true*y_sig + (1 - eta_true)*y_back


hs <- hist(obs, breaks = 30, probability = TRUE)

basic_plt<- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                              aes(y = after_stat(density)),
                                                            fill = 'steelblue', col = 'black',
                                                            breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = gb(xs), color = 'Proposal Background'),
            lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black')) + 
  geom_vline(xintercept = c(M_lower, M_upper),
             lwd = 1.3, lty = 2, color = 'cyan')

basic_plt
#-----------------------------------------------------------------

r_vec <- tau_vec <- aic_vec <- c()
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))
# obs_Gb <- sapply(obs, Gb)
for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs, T_basis[[j]])
  r_vec <- c(r_vec, cor(S1_vec, fun_vec))
  tau_vec <- c(tau_vec, mean(fun_vec))
  aic_vec <- c(aic_vec, mean(fun_vec)^2 - 2/n)
}

ord <- order(abs(tau_vec), decreasing = TRUE)

df <- data.frame(tau = tau_vec, r = r_vec, aic = aic_vec, 
                 order = ord, id = 1:length(T_basis))

df_ordered <- df[df$order,]
df_ordered$criteria <- abs(df_ordered$r)/sqrt(1-df_ordered$r^2)<qt(0.975, df = n)/sqrt(n)
# df_ordered$criteria <- sign(df_ordered$tau)*df_ordered$r<2/sqrt(n)
df_ordered$aic_cumsum <- cumsum(df_ordered$aic*df_ordered$criteria)

df_ordered

(aic_acc <- df_ordered$aic[df_ordered$criteria])
(fn_index <- df_ordered$id[df_ordered$criteria])

plot(y = aic_acc,
     x = seq_along(fn_index),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'AIC', xlab = 'Basis indices')
axis(1, at = seq_along(fn_index), labels = fn_index)

# final_tau <- c(0,0,tau_vec[3])

# final_tau <- c(rep(0,4),tau_vec[5]) # when bases on [l,u] is used

final_tau <- c(tau_vec[1]) # when bases on [l,u] with tuned gb is used

res <- sapply(xs, function(t) mod(t, theta = theta, tau = final_tau))

fitted_plot_plugins <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                                         aes(y = after_stat(density)),
                                                                       fill = 'steelblue', col = 'black',
                                                                       breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = gb(xs), color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[1,]), color = 'Fitted Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[2,]), color = 'Fitted Background from Mixture'),
            linetype = 'dashed', lwd = 1.3) + 
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black',
                                'Fitted Mixture' = 'cyan',
                                'Fitted Background from Mixture' = 'skyblue'))

fitted_plot_plugins

theta/norm_S

# testing for signal:

(se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n))
(t_stat <- theta/se_theta)
t_stat>qnorm(0.95)

