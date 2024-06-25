rm(list=ls())

source('functions_beta.R')
library(ggplot2)

data <- make_data(rate_back = rate_back,
                  a_sig = a_sig, b_sig = b_sig,
                  bkg_prop = 1-sig_prop,
                  seed = 12345)

xs <- seq(0,1,0.01)
ys <- dbeta(xs, shape1 = a_sig, shape2 = b_sig)
yb <- dtrunc(xs,spec="exp",rate=2,a=0,b=1)
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

r_vec <- NULL
mask <- NULL

for(i in 1:length(T_basis))
{
  Ti <- T_basis[[i]]
  Ti_vals <- sapply(mix, Ti, a = a_sig, b = b_sig)/norm_vec[i]
  S1_vals <- sapply(mix, 
                    function(t)
                      {
                      S(t, a = a_sig, b = b_sig)$s1
                    })
  r_val <- cor(Ti_vals, S1_vals)
  r_vec <- c(r_vec,
             abs(r_val))
  mask <- c(mask,
            abs(r_val)<=2/sqrt(2e3))
}

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
aic_terms <- tau1^2 - 2/2e3
bic_terms <- tau1^2 - log(2e3)/2e3
(AIC1_1 <- sum(aic_terms))
(BIC1_1 <- sum(bic_terms))
mask_1 <- mask[1:length(tau1)]
(AIC2_1 <- 0+sum(aic_terms[mask_1]))
(BIC2_1 <- 0+sum(bic_terms[mask_1]))


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
  geom_text(aes(label = paste0('AIC1: ', round(AIC1_1, 4),
                               '; BIC1: ',round(BIC1_1,4))), x=0.75, y = 2,
            size = 4) + 
  geom_text(aes(label = paste0('AIC2: ', round(AIC2_1, 4),
                               '; BIC2: ',round(BIC2_1,4))), x=0.75, y = 1.5,
            size = 4) + 
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
aic_terms <- tau2^2 - 2/2e3
bic_terms <- tau2^2 - log(2e3)/2e3
(AIC1_2 <- sum(aic_terms))
(BIC1_2 <- sum(bic_terms))
mask_2 <- mask[1:length(tau2)]
(AIC2_2 <- 0+sum(aic_terms[mask_2]))
(BIC2_2 <- 0+sum(bic_terms[mask_2]))


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
  geom_text(aes(label = paste0('AIC1: ', round(AIC1_2, 4),
                               '; BIC1: ',round(BIC1_2,4))), x=0.75, y = 2,
            size = 4)+ 
  geom_text(aes(label = paste0('AIC2: ', round(AIC2_2, 4),
                               '; BIC2: ',round(BIC2_2,4))), x=0.75, y = 1.5,
            size = 4)+ 
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
aic_terms <- tau3^2 - 2/2e3
bic_terms <- tau3^2 - log(2e3)/2e3
(AIC1_3 <- sum(aic_terms))
(BIC1_3 <- sum(bic_terms))
mask_3 <- mask[1:length(tau3)]
(AIC2_3 <- 0+sum(aic_terms[mask_3]))
(BIC2_3 <- 0+sum(bic_terms[mask_3]))



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
  geom_text(aes(label = paste0('AIC1: ', round(AIC1_3, 4),
                               '; BIC1: ',round(BIC1_3,4))), x=0.75, y = 2,
            size = 4) + 
  geom_text(aes(label = paste0('AIC2: ', round(AIC2_3, 4),
                               '; BIC2: ',round(BIC2_3,4))), x=0.75, y = 1.5,
            size = 4) + 
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
aic_terms <- tau4^2 - 2/2e3
bic_terms <- tau4^2 - log(2e3)/2e3
(AIC1_4 <- sum(aic_terms))
(BIC1_4 <- sum(bic_terms))
mask_4 <- mask[1:length(tau4)]
(AIC2_4 <- 0+sum(aic_terms[mask_4]))
(BIC2_4 <- 0+sum(bic_terms[mask_4]))


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
  geom_text(aes(label = paste0('AIC1: ', round(AIC1_4, 4),
                               '; BIC1: ',round(BIC1_4,4))), x=0.75, y = 2,
            size = 4) + 
  geom_text(aes(label = paste0('AIC2: ', round(AIC2_4, 4),
                               '; BIC2: ',round(BIC2_4,4))), x=0.75, y = 1.5,
            size = 4) + 
  ggtitle(bquote(T[1] + T[2] + T[3] + T[4])) +
  theme(plot.title = element_text(hjust = 0.5))
plt4

#-------------------------------------------------------------------------------
library(ggpubr)
final_plot <- ggarrange(plt1, plt2, plt3, plt4, nrow = 2, ncol = 2,
                                common.legend = TRUE, legend = 'bottom')
msg <- bquote(paste(
  .(1-sig_prop), '*', Exp[Trunc],'(',.(rate_back),')','+',
  .(sig_prop), '*',beta,'(',.(a_sig),',',.(b_sig),')' ))

annotate_figure(final_plot, top = text_grob(msg, size = 18, face = 'bold',
                                            color = 'red'))
