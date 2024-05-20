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
#---------------------------------------------------------------
t <- Sys.time()
res_MLE <- constrOptim(theta = rep(0.01,6),f = neg_loglikelihood,
                       data = obs,
                       grad = score, 
                       ui = rbind(c(1,0,0,0,0,0), c(-1,0,0,0,0,0)),
                       ci = c(0,-norm_S))
(time_taken <- Sys.time() - t)

tau_MLE <- res_MLE$par[-1]
theta_MLE <- res_MLE$par[1]

ord_MLE <- order(abs(tau_MLE), decreasing = TRUE)

df <- data.frame(tau = tau_MLE, order = ord_MLE, 
                 id = 1:length(T_basis))
df_ordered <- df[df$order,]
df_ordered
# df_ordered$criteria <- abs(df_ordered_MLE$r)<2/sqrt(n)
aic_vec <- bic_vec <- c()
for(j in 1:length(tau_MLE))
{
  temp_tau <- rep(0, length(tau_MLE))
  temp_tau[df_ordered$id[1:j]] <- df_ordered$tau[1:j]
  temp_params <- c(theta_MLE, temp_tau)
  aic_vec <- c(aic_vec,
               (-2)*neg_loglikelihood(temp_params, data = obs) - 2*j)
  bic_vec <- c(bic_vec,
               (-2)*neg_loglikelihood(temp_params, data = obs) - log(n)*j)
}
df_ordered$aic <- aic_vec - min(aic_vec)
df_ordered$bic <- bic_vec - min(bic_vec)

par(mfrow = c(1,2))
plot(y = df_ordered$aic,
     x = seq_along(df_ordered$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'AIC', xlab = 'Basis indices')
axis(1, at = seq_along(df_ordered$id), labels = df_ordered$id)

plot(y = df_ordered$bic,
     x = seq_along(df_ordered$id),type = 'o', lwd = 3, xaxt = "none",
     ylab = 'BIC', xlab = 'Basis indices')
axis(1, at = seq_along(df_ordered$id), labels = df_ordered$id)

final_tau_MLE <- tau_MLE[c(1,3)]

res_MLE_final <- sapply(xs, function(t) mod(t, theta = theta_MLE, tau = final_tau_MLE))

fitted_plot_MLE <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                                     aes(y = after_stat(density)),
                                                                   fill = 'steelblue', col = 'black',
                                                                   breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_real, color = 'Actual Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_back, color = 'Background'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = gb(xs), color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res_MLE_final[1,]), color = 'Fitted Mixture-MLE'),
            linetype = 'dashed', lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res_MLE_final[2,]), color = 'Fitted Background from Mixture-MLE'),
            linetype = 'dashed', lwd = 1.3) +
  scale_color_manual(values = c('Background' = 'orange',
                                'Proposal Background' = 'red',
                                'Actual Mixture' = 'black',
                                'Fitted Mixture-MLE' = 'cyan',
                                'Fitted Background from Mixture-MLE' = 'skyblue'))

fitted_plot_MLE

