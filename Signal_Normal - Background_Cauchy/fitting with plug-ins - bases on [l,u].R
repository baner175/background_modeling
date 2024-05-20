library(ggplot2)
source('bases on [l,u].R')

set.seed(12345)

n <- 5e3

sig <- rtrunc(n, spec = 'norm',
              mean = mean_sig, sd = sd_sig,
              a = l, b = u)
back <- rtrunc(n, spec = 'cauchy',
               location = mean_back, scale = sd_back,
               a = l, b = u)
u_mask <- runif(n)
obs <- ifelse(u_mask<1-eta_true,back,sig)

xs <- seq(l,u,0.01)
y_sig <- dtrunc(xs, spec = 'norm',
                mean = mean_sig, sd = sd_sig,
                a = l, b = u)
y_back <- dtrunc(xs, spec = 'cauchy', 
                 location = mean_back, scale = sd_back,
                 a = l, b = u)
y_real <- eta_true*y_sig + (1 - eta_true)*y_back


hs <- hist(obs, breaks = 20, probability = TRUE)

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

# checking conditions----------------------------------------------
# checking delta<0:
delta <- integrate(function(t) S1(t, 
                         mean = mean_sig, sd = sd_sig)* fb(t), l, u)$value

delta

#--------------------------------------------------------

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

df_ordered$criteria <- abs(df_ordered$r)<2/sqrt(n)
# df_ordered$criteria <- sign(df_ordered$tau)*df_ordered$r<2/sqrt(n)
df_ordered$aic_cumsum <- cumsum(df_ordered$aic*df_ordered$criteria)

df_ordered

aic_conv <- bic_conv <- c()
for(j in 1:length(tau_vec))
{
  temp_tau <- rep(0, length(tau_vec))
  temp_tau[df_ordered$id[1:j]] <- df_ordered$tau[1:j]
  temp_params <- c(theta, temp_tau)
  aic_conv <- c(aic_conv,
                (-2)*neg_loglikelihood(temp_params, data = obs) - 2*j)
  bic_conv <- c(bic_conv,
                (-2)*neg_loglikelihood(temp_params, data = obs) - log(n)*j)
}
df_ordered$aic_conv <- aic_conv
df_ordered$bic_conv <- bic_conv

par(mfrow = c(1,3))
plot(y = df_ordered$aic_cumsum,
     x = seq_along(df_ordered$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'AIC', xlab = 'Basis indices')
axis(1, at = seq_along(df_ordered$id), labels = df_ordered$id)

plot(y = df_ordered$aic_conv,
     x = seq_along(df_ordered$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'AIC - conventional', xlab = 'Basis indices')
axis(1, at = seq_along(df_ordered$id), labels = df_ordered$id)

plot(y = df_ordered$bic_conv,
     x = seq_along(df_ordered$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'BIC - conventional', xlab = 'Basis indices')
axis(1, at = seq_along(df_ordered$id), labels = df_ordered$id)

final_tau <- tau_vec[c(1,2,3)]

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
