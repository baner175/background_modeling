# Estimated eta:  0.01837047
# Test for eta: null hypothesis (eta = 0) rejected with p-value 0.01040974

library(ggplot2)

source('bases on [l,u] with tuned gb.R')
dat <- read.table('Data_ex1.txt', header = TRUE)
obs <- log(dat$x)
n <- length(obs)

xs <- seq(l, u, length.out = 1e3)
y_sig <- sapply(xs, fs)
y_gb <- sapply(xs, gb)

kde <- kdensity::kdensity(obs)

hs <- hist(obs, probability = TRUE, breaks = 50)
basic_plt <- ggplot(mapping = aes(x = obs)) + 
  geom_histogram(aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(aes(x = xs, y = y_gb), col = 'red',
            lwd = 1.2, alpha = 0.5) +
  geom_line(aes(x = xs, y = kde(xs)), col = 'yellow',
            lwd = 1.2, alpha = 0.5) + 
  geom_vline(xintercept = c(log(mean_sig), M_lower, M_upper), 
             col = c('black', 'blue', 'blue'),
             lty = 2, lwd = 1.2, alpha = 0.5)
basic_plt

#-----------------------------------------------------------------

r_vec <- tau_vec <- se_vec <- c()
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))

se_theta <- sqrt((mean(S1_vec^2) - theta^2)/n)

# testing \eta = 0:
theta/norm_S

# testing for signal:
(t_stat_theta <- theta/se_theta)
pnorm(t_stat_theta, lower.tail = FALSE)
t_stat_theta>qnorm(0.95)

(theta <- theta*(t_stat_theta>qnorm(0.95)))

for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs, T_basis[[j]])
  r_vec <- c(r_vec, cor(S1_vec, fun_vec))
  tau_vec <- c(tau_vec, mean(fun_vec))
  se_vec <- c(se_vec, sqrt((mean(fun_vec^2) - mean(fun_vec)^2)/n))
}


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
            tau = t) - 2*(sum(t!=0) + (theta!=0))
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

fitted_plot_plugins <- ggplot(mapping = aes(x = obs)) + geom_histogram(mapping = 
                                                                         aes(y = after_stat(density)),
                                                                       fill = 'steelblue', col = 'black',
                                                                       breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = gb(xs), color = 'Proposal Background'),
            lwd = 1.3, alpha = 0.8) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[1,]), color = 'Fitted Mixture'),
            linetype = 'dashed', lwd = 1.3, alpha = 0.8) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[2,]), color = 'Fitted Background from Mixture'),
            linetype = 'dashed', lwd = 1.3, alpha = 0.8) + 
  geom_line(aes(x = xs, y = kde(xs), color = 'KDE for mixture'),
            lwd = 1.2, alpha = 0.8) +
  scale_color_manual(values = c('Proposal Background' = 'red',
                                'Fitted Mixture' = 'cyan',
                                'Fitted Background from Mixture' = 'skyblue',
                                'KDE for mixture' = 'yellow'))

fitted_plot_plugins

theta/norm_S

se_theta/norm_S


# The plot containing both unknown and estimated densities
phi_bkg <- 1.4
eta_true <- 0.02

fb <- function(x)
{
  return(dtrunc(exp(x), spec = 'pareto', 
                shape = phi_bkg, scale = real_l,
                a = real_l, b = real_u)*exp(x))
}

f_mix <- function(x)
{
  return((1-eta_true)*fb(x)+eta_true*fs(x))
}
y_bkg <- sapply(xs, fb)
y_mix <- sapply(xs, f_mix)

My_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(mapping = aes(x = obs)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'grey', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_gb),
            lwd = 1.2, col = 'orange', linetype = 6) + 
  geom_line(mapping = aes(x = xs, y = y_bkg),
            lwd = 1.2, col = 'red', linetype = 2) + 
  geom_line(mapping = aes(x = xs, y = y_mix),
            lwd = 1.2, col = 'blue', linetype = 1) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[1,])),
            col = 'cyan',linetype = 4, lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(res[2,])), 
            col = 'brown',linetype = 5, lwd = 1.2) + 
  xlab('log(y)') + ylab('Density') +
  theme_bw() +
  My_Theme


#Checking delta:

integrate(function(t) S1(t)*fb(t), l, u)
