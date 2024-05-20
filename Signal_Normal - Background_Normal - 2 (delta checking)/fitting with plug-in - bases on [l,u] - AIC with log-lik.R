source('bases on [l,u].R')

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
                                'Actual Mixture' = 'black'))

basic_plt
#-----------------------------------------------------------------

r_vec <- tau_vec <- aic_vec <- c()
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))
# obs_Gb <- sapply(obs, Gb)
for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs, T_basis[[j]])
  tau_vec <- c(tau_vec, mean(fun_vec))
}

ord <- order(abs(tau_vec), decreasing = TRUE)

df <- data.frame(tau = tau_vec, id = 1:length(T_basis))


aic_vec <- bic_vec <- c()

for(j in 1:length(tau_vec))
{
  temp_tau <- df$tau[1:j]
  aic_vec <- c(aic_vec,
               2*log_lik(data = obs, theta = theta, 
                         tau = temp_tau) - 2*j)
  bic_vec <- c(bic_vec,
               2*log_lik(data = obs, theta = theta, 
                         tau = temp_tau) - log(n)*j)
}

df$aic <- aic_vec - min(aic_vec)
df$bic <- bic_vec - min(bic_vec)

par(mfrow = c(1,2))
plot(y = df$aic,
     x = seq_along(df$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'AIC', xlab = 'Basis indices')
axis(1, at = seq_along(df$id), labels = df$id)

plot(y = df$bic,
     x = seq_along(df$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'BIC', xlab = 'Basis indices')
axis(1, at = seq_along(df$id), labels = df$id)



# df_ordered$criteria <- sign(df_ordered$tau)*df_ordered$r<2/sqrt(n)

final_tau <- tau_vec[1:2]

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
