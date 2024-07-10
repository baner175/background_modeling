# ESTIMATED eta: 0.02319415
# Test for eta: null hypothesis (eta = 0) rejected with p-value 1.727512e-05


source('bases on [l,u] - gb with 2 bumps.R')

library(ggplot2)

xs <- seq(l,u,0.01)
y_sig <- sapply(xs, fs)
y_kde <- sapply(xs, kde)
y_gb <- sapply(xs, gb)

obs_new <- read.csv('full_data.csv', header = TRUE)$x
hs <- hist(obs_new, breaks = 50, probability = TRUE)

basic_plt<- ggplot(mapping = aes(x = obs_new)) + 
  ylim(0, .5) +
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_kde, color = 'KDE for Mixture'),
            lwd = 1.3) +
  geom_line(mapping = aes(x = xs, y = y_gb, color = 'Proposal Background'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = y_sig, color = 'True Signal'),
            lwd = 1.3) + 
  scale_color_manual(values = c('KDE for Mixture' = 'black',
                                'Proposal Background' = 'red',
                                'True Signal' = 'blue')) + 
  geom_vline(xintercept = c(M_lower, M_upper),
             lwd = 1.3, lty = 2, color = 'grey')

basic_plt
#-----------------------------------------------------------------

r_vec <- tau_vec <- se_vec <- c()
S1_vec <- sapply(obs_new, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))
se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)

for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs_new, T_basis[[j]])
  r_vec <- c(r_vec, cor(S1_vec, fun_vec))
  tau_vec <- c(tau_vec, mean(fun_vec))
  se_vec <- c(se_vec, sqrt((mean(fun_vec^2) - mean(fun_vec)^2)/n))
}

# testing \eta = 0:
theta/norm_S

se_theta/norm_S

# testing for signal:
(t_stat_theta <- theta/se_theta)
t_stat_theta>qnorm(0.95)

pnorm(t_stat_theta, lower.tail = FALSE)

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
  2*log_lik(data = obs_new, theta = theta, 
            tau = t) - 2*(sum(t!=0) + (theta != 0))
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

(final_tau <- taus_for_aic[which.max(df$aic),])

res <- sapply(xs, function(t) mod(t, theta = theta, tau = final_tau))


#----------------------------------------------------------------

mean_back <- 0.5; sd_back <- 2.5
eta_true <- 0.03
f_back <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                             mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*f_back(t)

y_mix <- sapply(xs, f_mix)
y_back <- sapply(xs, f_back)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(mapping = aes(x = obs)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 breaks = hs$breaks, col = 'black',
                 fill = 'gray') + 
  geom_line(mapping = aes(x = xs, y = y_mix), col = 'blue',
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs, y = y_back), col = 'red',
            lwd = 1.3, linetype = 2) + 
  geom_line(mapping = aes(x = xs, y = y_gb), col = 'orange',
            lwd = 1.3, linetype = 6) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[1,])),
            col = 'cyan',linetype = 4, lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[2,])), 
            col = 'brown',linetype = 5, lwd = 1.2) + 
  xlab('y') + ylab('Density') +
  theme_bw() +
  My_Theme


mean(S1_vec)/norm_S

#Checking delta:

integrate(function(t) S1(t)*f_back(t), l, u)
