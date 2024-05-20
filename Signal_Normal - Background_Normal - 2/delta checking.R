source('bases on [l,u].R')

fb <- function(x) dtrunc(x, spec = 'norm', 
                         mean = mean_back, sd = sd_back,
                         a = l, b = u)

f <- function(x) eta_true*fs(x) + (1-eta_true)*fb(x)

library(ggplot2)

set.seed(12345)

# checking if delta < 0:
integrate(function(t) S1(t, mean = mean_sig,
                         sd = sd_sig)*dtrunc(t, mean = mean_back,
                                             sd = sd_back, spec = 'norm',
                                             a = l, b = u),
          lower = l, upper = u)$value

n <- 5e3

back <- rtrunc(n, spec = 'norm',
               mean = mean_back, sd = sd_back,
               a = l, b = u)
sig <- rtrunc(n, spec = 'norm',
              mean = mean_sig, sd = sd_sig,
              a = l, b = u)
u_mask <- rbinom(n, size = 1, prob = eta_true)

obs <- ifelse(u_mask, sig, back)

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
#-------------------------------------------------------------------------------

tau_true <- c()

for(j in 1:length(T_basis))
{
  real_inner <- integrate(function(t) T_basis[[j]](t)*f(t), l, u)$value
  tau_true <- c(tau_true, real_inner)
}
tau_true
(theta_true <- integrate(function(t) S1(t)*f(t), l, u)$value)
(delta_true <- integrate(function(t) S1(t)*fb(t), l, u)$value)

sum_tau_true <- sum(tau_true^2)
(eta_tilde_true <- theta_true/norm_S)

LHS_num_true <- sum_tau_true + norm_S^2 - norm_S*theta_true
LHS_den_true <- sqrt(sum_tau_true + theta_true^2 + norm_S^2 - 2*norm_S*theta_true)
LHS_true <- LHS_num_true/LHS_den_true

RHS_num_true <- sum_tau_true/(1-eta_tilde_true) + norm_S^2
RHS_den_true <- sqrt(sum_tau_true/((1-eta_tilde_true)^2) + norm_S^2 )
RHS_true <- RHS_num_true/RHS_den_true

LHS_true>RHS_true

round(LHS_true, 10) == round(RHS_true, 10)

(norm_fb_perp_bar_true <- sqrt(sum(tau_true^2)/((1-eta_tilde_true)^2)))
(norm_fb_perp_true <- sqrt(sum(tau_true^2)/((1-eta_true)^2)))
(norm_f_perp_true <- sqrt(sum(tau_true^2)))


plot(x = seq(-0.001,norm_fb_perp_bar_true+0.0001,length.out=10),
     y = seq(delta_true-0.0001, norm_S+1, length.out = 10),
     type = 'n',
     xlab = '', ylab = '')
abline(h = 0, lwd = 2, col = 'red')
lines(x = c(0, 0), y = c(0, norm_S), col = 'blue', lwd = 2)
lines(x = c(norm_fb_perp_true, norm_fb_perp_true), y = c(0, delta_true),
      col = 'brown', lwd = 2)
lines(x = c(0, norm_fb_perp_true),
      y = c(norm_S, delta_true),
      lwd = 2, col = 'black')
lines(x = c(norm_f_perp_true, norm_f_perp_true),
      y = c(0, theta_true))
lines(x = c(0,norm_f_perp_true),
      y = c(norm_S, 0), 
      col = 'green')
abline(v = norm_fb_perp_bar_true, lty = 2)
lines(x = c(0,norm_fb_perp_true),
      y = c(norm_S, 0), 
      col = 'green')



#-------------------------------------------------------------------------------

tau_vec <- c()
S1_vec <- sapply(obs, function(t) S1(t, mean = mean_sig,
                                     sd = sd_sig))
(theta <- mean(S1_vec))
# obs_Gb <- sapply(obs, Gb)
for (j in 1:length(T_basis)) {
  fun_vec <- sapply(obs, T_basis[[j]])
  tau_vec <- c(tau_vec, mean(fun_vec))
}

sum_tau <- sum(tau_vec^2)
eta_tilde <- theta/norm_S

LHS_num <- sum_tau + norm_S^2 - norm_S*theta
LHS_den <- sqrt(sum_tau + theta^2 + norm_S^2 - 2*norm_S*theta)
LHS <- LHS_num/LHS_den

RHS_num <- sum_tau/(1-eta_tilde) + norm_S^2
RHS_den <- sqrt( sum_tau/((1-eta_tilde)^2) + norm_S^2 )
RHS <- RHS_num/RHS_den

LHS

RHS

round(LHS, 10) == round(RHS, 10)
#-------------------------------------------------------------------------------

(delta_true <- integrate(function(t) S1(t)*fb(t), l, u)$value)

plot(x = seq(-0.001,sqrt(sum(tau_vec^2)/((1-eta_tilde)^2))+0.0001,length.out=10),
     y = seq(delta_true+0.00001, norm_S+1, length.out = 10),
     type = 'n',
     xlab = '', ylab = '')
abline(h = 0, lwd = 2, col = 'red')
lines(c(sqrt(sum(tau_vec^2)/((1-eta_tilde)^2)),sqrt(sum(tau_vec^2)/((1-eta_tilde)^2))),
      c(0,delta_true))

lines(c(0, 0),
      c(0,norm_S))
lines()