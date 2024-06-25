rm(list=ls())
library(ggplot2)
library(truncdist)
library(VGAM)
rate_sig <- 2
sig_prop <- 0.05
back_scale <- 2; back_shape <- 2

source('functions_exp.R')

set.seed(12345)
bkg<- rtrunc(2e3,spec = 'pareto',
             scale = back_scale, 
             shape = back_shape, a = 2,b=3)
sig <- rtrunc(2e3, spec = 'exp', 
              rate = rate_sig,
              a = 2, b = 3)
u_mask <- runif(2e3)
mix <- ifelse(u_mask<sig_prop, sig, bkg)

xs <- seq(2,3,0.01)
ys <- dtrunc(xs, spec = 'exp', 
             rate = rate_sig,
             a = 2, b = 3)
yb <- dtrunc(xs,spec = 'pareto',
             scale = back_scale, 
             shape = back_shape, a = 2,b=3)
y_mixed <- ys*sig_prop + yb*(1-sig_prop)

hs <- hist(mix, probability = TRUE, breaks = 20)
ggplot(mapping = aes(x = mix)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue',
                 col = 'black', breaks = hs$breaks) +
  geom_line(mapping = aes(x = xs,
                          y = ys, color = 'sig'), 
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs,
                          y = yb, color = 'bkg'), 
            lwd = 1.3) + 
  geom_line(mapping = aes(x= xs,
                          y = 1, color = 'proposal bkg'),
            lwd = 1.3) + 
  geom_line(mapping = aes(x = xs,
                          y = y_mixed, color = 'mixture'),
            lwd = 1.3) + 
  scale_color_manual(values = c('bkg' = 'orange',
                                'sig' = 'skyblue',
                                'proposal bkg' = 'red',
                                'mixture' = 'black'))

# r_vec <- NULL
# mask <- NULL
# 
# for(i in 1:length(T_basis))
# {
#   Ti <- T_basis[[i]]
#   Ti_vals <- sapply(mix, Ti, rate = rate_sig)/norm_vec[i]
#   S1_vals <- sapply(mix, 
#                     function(t)
#                     {
#                       S(t, rate = rate_sig)$s1
#                     })
#   r_val <- cor(Ti_vals, S1_vals)
#   r_vec <- c(r_vec,
#              abs(r_val))
#   mask <- c(mask,
#             abs(r_val)<=2/sqrt(2e3))
# }

#------------------------------------------------------------------------
res1 <- constrOptim(theta = rep(0.1,2),f = neg_loglikelihood,
                    data = mix, rate_sig = rate_sig,
                    grad = score, 
                    ui = rbind(c(1,0),c(-1,0)),
                    ci = c(0,-1))

back_mle_1 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  rate_sig = rate_sig, p = 0.01)

(eta <- res1$par[1])
(beta <- res1$par[-1])
(beta_back <- back_mle_1$estimate)
ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, rate = rate_sig)
})

bkml_1 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, rate = rate_sig)
})



plt1 <- ggplot(mapping = aes(x = mix)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_mixed, color = 'True Mixture'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys1[1,]), 
                          color = 'Estimated Mixture'),
            lwd = 1.2, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys1[2,]), 
                          color = 'Estimated Background from Mixture'),
            lwd = 1.2, linetype = 'dotted') + 
  geom_line(mapping = aes(x = xs, y = yb, 
                          color = 'True Background'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = bkml_1, 
                          color = 'Estmated Background form Background only'),
            lwd = 1.2, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = 1, 
                          color = 'Proposed Background'),
            lwd = 1.2) +
  scale_color_manual(values = c('True Mixture' = 'black',
                                'Estimated Mixture' = 'blue',
                                'Estimated Background from Mixture' = 'green',
                                'True Background' = 'brown',
                                'Estmated Background form Background only' = 'orange',
                                'Proposed Background' = 'red')) +
  ggtitle(bquote(T[1])) + 
  theme(plot.title = element_text(hjust = 0.5))

plt1
#-------------------------------------------------------------------------------

res2 <- constrOptim(theta = rep(0.001,3),f = neg_loglikelihood,
                    data = mix, rate_sig = rate_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))

back_mle_2 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  rate_sig = rate_sig, p = rep(0.01,2))

(eta <- res2$par[1])
(beta <- res2$par[-1])
(beta_back <- back_mle_2$estimate)
ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, rate = rate_sig)
})

bkml_2 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, rate = rate_sig)
})



plt2 <- ggplot(mapping = aes(x = mix)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_mixed, color = 'True Mixture'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys2[1,]), 
                          color = 'Estimated Mixture'),
            lwd = 1.2, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys2[2,]), 
                          color = 'Estimated Background from Mixture'),
            lwd = 1.2, linetype = 'dotted') + 
  geom_line(mapping = aes(x = xs, y = yb, 
                          color = 'True Background'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = bkml_2, 
                          color = 'Estmated Background form Background only'),
            lwd = 1.2, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = 1, 
                          color = 'Proposed Background'),
            lwd = 1.2) +
  scale_color_manual(values = c('True Mixture' = 'black',
                                'Estimated Mixture' = 'blue',
                                'Estimated Background from Mixture' = 'green',
                                'True Background' = 'brown',
                                'Estmated Background form Background only' = 'orange',
                                'Proposed Background' = 'red')) +
  ggtitle(bquote(T[1] + T[2])) + 
  theme(plot.title = element_text(hjust = 0.5))
plt2
#-------------------------------------------------------------------------------

res3 <- constrOptim(theta = rep(0.001,4),f = neg_loglikelihood,
                    data = mix, rate_sig = rate_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))

back_mle_3 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  rate_sig = rate_sig, p = rep(0.01,3))

(eta <- res3$par[1])
(beta <- res3$par[-1])
(beta_back <- back_mle_3$estimate)
ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, rate = rate_sig)
})

bkml_3 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, rate = rate_sig)
})



plt3 <- ggplot(mapping = aes(x = mix)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_mixed, color = 'True Mixture'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys3[1,]), 
                          color = 'Estimated Mixture'),
            lwd = 1.2, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys3[2,]), 
                          color = 'Estimated Background from Mixture'),
            lwd = 1.2, linetype = 'dotted') + 
  geom_line(mapping = aes(x = xs, y = yb, 
                          color = 'True Background'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = bkml_3, 
                          color = 'Estmated Background form Background only'),
            lwd = 1.2, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = 1, 
                          color = 'Proposed Background'),
            lwd = 1.2) +
  scale_color_manual(values = c('True Mixture' = 'black',
                                'Estimated Mixture' = 'blue',
                                'Estimated Background from Mixture' = 'green',
                                'True Background' = 'brown',
                                'Estmated Background form Background only' = 'orange',
                                'Proposed Background' = 'red')) + 
  ggtitle(bquote(T[1] + T[2] + T[3])) +
  theme(plot.title = element_text(hjust = 0.5))
plt3
#-------------------------------------------------------------------------------

res4 <- constrOptim(theta = rep(0.001,5),f = neg_loglikelihood,
                    data = mix, rate_sig = rate_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0,0),c(-1,0,0,0,0)),
                    ci = c(0,-1))

back_mle_4 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  rate_sig = rate_sig, p = rep(0.01,4))

(eta <- res4$par[1])
(beta <- res4$par[-1])
(beta_back <- back_mle_4$estimate)
ys4 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, rate = rate_sig)
})

bkml_4 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, rate = rate_sig)
})


plt4 <- ggplot(mapping = aes(x = mix)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_mixed, color = 'True Mixture'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys4[1,]), 
                          color = 'Estimated Mixture'),
            lwd = 1.2, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys4[2,]), 
                          color = 'Estimated Background from Mixture'),
            lwd = 1.2, linetype = 'dotted') + 
  geom_line(mapping = aes(x = xs, y = yb, 
                          color = 'True Background'),
            lwd = 1.2) +
  geom_line(mapping = aes(x = xs, y = bkml_3, 
                          color = 'Estmated Background form Background only'),
            lwd = 1.2, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = 1, 
                          color = 'Proposed Background'),
            lwd = 1.2) +
  scale_color_manual(values = c('True Mixture' = 'black',
                                'Estimated Mixture' = 'blue',
                                'Estimated Background from Mixture' = 'green',
                                'True Background' = 'brown',
                                'Estmated Background form Background only' = 'orange',
                                'Proposed Background' = 'red')) + 
  ggtitle(bquote(T[1] + T[2] + T[3] + T[4])) +
  theme(plot.title = element_text(hjust = 0.5))
plt4

#-------------------------------------------------------------------------------
library(ggpubr)
final_plot <- ggarrange(plt1, plt2, plt3, nrow = 1, ncol = 3,
                        common.legend = TRUE, legend = 'bottom')
msg <- bquote(paste(
  .(sig_prop), '*', Exp[Trunc],'(',.(rate_sig),')','+',
  .(1-sig_prop), '*',Pareto, '(',.(back_scale),',',.(back_shape),')'))

annotate_figure(final_plot, top = text_grob(msg, size = 18, face = 'bold',
                                            color = 'red'))
