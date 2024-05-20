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

r_vec <- tau_vec <- se_vec <- c()
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))
se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)

for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs, T_basis[[j]])
  r_vec <- c(r_vec, cor(S1_vec, fun_vec))
  tau_vec <- c(tau_vec, mean(fun_vec))
  se_vec <- c(se_vec, sqrt((mean(fun_vec^2) - mean(fun_vec)^2)/n))
}

# testing \eta = 0:
theta/norm_S

# testing for signal:
(t_stat_theta <- theta/se_theta)
t_stat_theta>qnorm(0.95)

(theta <- theta*(t_stat_theta>qnorm(0.95)))

corr_test_taus <- abs(r_vec)/sqrt(1-r_vec^2)<qt(0.975, df = n)/sqrt(n)

mask_for_funs <- c(corr_test_taus)

coefs_to_use <- c(tau_vec)*mask_for_funs


df <- data.frame(fun_names = paste0('T',1:5),
                 coefs = coefs_to_use,
                 T_stat = c(coefs_to_use/se_vec))
df <- df[order(abs(df$T_stat), decreasing = TRUE),]

df

fun_indices <- as.numeric(sapply(strsplit(df$fun_names, 'T'), function(t) t[2]))
fun_indices <- fun_indices[df$coefs != 0]

taus_for_aic <- c()
for(i in 1:length(fun_indices))
{
  vec <- rep(0, length(T_basis))
  vec[fun_indices[1:i]] <- 1
  taus_for_aic <- rbind(taus_for_aic, vec*tau_vec)
}

aic_vec <- apply(taus_for_aic, 1, function(t){
  2*log_lik(data = obs, theta = theta, 
            tau = t) - 2*sum(t!=0)
  })

aic_vec <- c(aic_vec, rep(NA, nrow(df) - length(aic_vec)))

df$aic <- aic_vec

custom_order <- df$fun_names

df$fun_names <- factor(df$fun_names,
                       levels = custom_order)

ggplot(df[complete.cases(df),], aes(x = fun_names, y = aic)) + 
  geom_point(size = 3, col = 'blue', alpha = 0.5) + 
  geom_line(mapping = aes(y = aic,
                          x = 1:length(aic)), lwd = 1.2,
            alpha = 0.5)
  
final_tau <- taus_for_aic[which.max(df$aic),]

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
