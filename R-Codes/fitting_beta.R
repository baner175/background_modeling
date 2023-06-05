rm(list=ls())

source('functions_beta.R')

data <- make_data(a_back = a_back, b_back = b_back,
                  a_sig = a_sig, b_sig = b_sig,
                  bkg_prop = 1-sig_prop,
                  seed = 12345)

xs <- seq(0,1,0.01)
ys <- dbeta(xs, shape1 = a_sig, shape2 = b_sig)
yb <- dbeta(xs, shape1 = a_back, shape2 = b_back)
y_mixed <- ys*sig_prop + yb*(1-sig_prop)

sig <- data$signal
bkg <- data$background
mix <- data$observed

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

#------------------------------------------------------------------------
res1 <- constrOptim(theta = rep(0.1,2),f = neg_loglikelihood,
                    data = mix, a_sig = a_sig, b_sig = b_sig,
                    grad = score, 
                    ui = rbind(c(1,0),c(-1,0)),
                    ci = c(0,-1))

back_mle_1 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  a_sig = a_sig, b_sig = b_sig, p = 0.01)

(eta <- res1$par[1])
(beta <- res1$par[-1])
(beta_back <- back_mle_1$estimate)
ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, a = a_sig, b = b_sig)
})

bkml_1 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, a = a_sig, b = b_sig)
})

tau1 <- res1$par[-1]*(1-res1$par[1])
(AIC1 <- sum(tau1^2) - (2*1/2e3))
(BIC1 <- sum(tau1^2) - (1*log(2e3)/2e3))

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
  geom_text(aes(label = paste0('AIC: ', round(AIC1, 4),
                               '; BIC: ',round(BIC1,4))), x=0.5, y = 5,
            size = 5) + 
  ggtitle(bquote(T[1])) + 
  theme(plot.title = element_text(hjust = 0.5))

plt1
#-------------------------------------------------------------------------------

res2 <- constrOptim(theta = rep(0.1,3),f = neg_loglikelihood,
                    data = mix, a_sig = a_sig, b_sig = b_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))

back_mle_2 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  a_sig = a_sig, b_sig = b_sig, p = rep(0.01,2))

(eta <- res2$par[1])
(beta <- res2$par[-1])
(beta_back <- back_mle_2$estimate)
ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, a = a_sig, b = b_sig)
})

bkml_2 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, a = a_sig, b = b_sig)
})

tau2 <- res2$par[-1]*(1-res2$par[1])
(AIC2 <- sum(tau2^2) - (2*2/2e3))
(BIC2 <- sum(tau2^2) - (2*log(2e3)/2e3))


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
  geom_text(aes(label = paste0('AIC: ', round(AIC2, 4),
                               '; BIC: ',round(BIC2,4))), x=0.5, y = 5,
            size = 5)+ 
  ggtitle(bquote(T[1] + T[2])) + 
  theme(plot.title = element_text(hjust = 0.5))
plt2
#-------------------------------------------------------------------------------

res3 <- constrOptim(theta = rep(0.1,4),f = neg_loglikelihood,
                    data = mix, a_sig = a_sig, b_sig = b_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))

back_mle_3 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  a_sig = a_sig, b_sig = b_sig, p = rep(0.01,3))

(eta <- res3$par[1])
(beta <- res3$par[-1])
(beta_back <- back_mle_3$estimate)
ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, a = a_sig, b = b_sig)
})

bkml_3 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, a = a_sig, b = b_sig)
})

tau3 <- res3$par[-1]*(1-res3$par[1])
(AIC3 <- sum(tau3^2) - (2*3/2e3))
(BIC3 <- sum(tau3^2) - (3*log(2e3)/2e3))


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
  ylim(-1,5) + 
  geom_text(aes(label = paste0('AIC: ', round(AIC3, 4),
                               '; BIC: ',round(BIC3,4))), x=0.5, y = 5,
            size = 5) + 
  ggtitle(bquote(T[1] + T[2] + T[3])) +
  theme(plot.title = element_text(hjust = 0.5))
plt3
#-------------------------------------------------------------------------------

res4 <- constrOptim(theta = rep(0.1,5),f = neg_loglikelihood,
                    data = mix, a_sig = a_sig, b_sig = b_sig,
                    grad = score, 
                    ui = rbind(c(1,0,0,0,0),c(-1,0,0,0,0)),
                    ci = c(0,-1))

back_mle_4 <- nlm(f = neg_loglikelihood_back_new, data = bkg,
                  a_sig = a_sig, b_sig = b_sig, p = rep(0.01,4))

(eta <- res4$par[1])
(beta <- res4$par[-1])
(beta_back <- back_mle_4$estimate)
ys4 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, a = a_sig, b = b_sig)
})

bkml_4 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, a = a_sig, b = b_sig)
})

tau4 <- res4$par[-1]*(1-res4$par[1])
(AIC4 <- sum(tau4^2) - (2*4/2e3))
(BIC4 <- sum(tau2^2) - (4*log(2e3)/2e3))


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
  ylim(-1,5) +
  geom_text(aes(label = paste0('AIC: ', round(AIC4, 4),
                               '; BIC: ',round(BIC4,4))), x=0.5, y = 5,
            size = 5) + 
  ggtitle(bquote(T[1] + T[2] + T[3] + T[4])) +
  theme(plot.title = element_text(hjust = 0.5))
plt4

#-------------------------------------------------------------------------------
final_plot <- ggpubr::ggarrange(plt1, plt2, plt3, plt4, nrow = 2, ncol = 2,
                                common.legend = TRUE, legend = 'bottom')
msg <- bquote(paste(
  .(1-sig_prop), '*', beta,'(',.(a_back),',',.(b_back),')','+',
  .(sig_prop), '*',beta,'(',.(a_sig),',',.(b_sig),')' ))

annotate_figure(final_plot, top = text_grob(msg, size = 18, face = 'bold',
                                            color = 'red'))
