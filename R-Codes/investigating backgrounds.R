rm(list = ls())

source('functions.R')
source('bkg_MLE functions.R')

library(patchwork)
library(ggplot2)

par(mfrow = c(1,3))


data <- make_data(mean_sig = 0.4, 
                  sd_sig = 0.1,
                  mean_back = 5,
                  sd_back = 5,
                  bkg_prop = 0.97,
                  seed = 12345)

data_back <- data$background
data <- data$observed

xs <- seq(0,1,0.01)
y_real <- sapply(xs, actual, mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97)


res1 <- constrOptim(theta = rep(0.1,2),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0),c(-1,0)),
                    ci = c(0,-1))

back_mle1_new <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = 0.01)

back_mle1_old <- nlm(f = neg_loglikelihood_back, data = data_back, 
                     p = 0.01)
back_mle1_old_full <- nlm(f = neg_loglikelihood_back, data = data, 
                      p = 0.01)
(eta <- res1$par[1])
(beta <- res1$par[-1])
(beta_back_new <- back_mle1_new$estimate)
(beta_back_old <- back_mle1_old$estimate)
(beta_back_old_full <- back_mle1_old_full$estimate)

ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml_new1 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back_new, mean = 0.4, sd = 0.1)
})

bkml_old1 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old)
})

bkml_old_full1 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old_full)
})

colors <- c('mixed data, new basis' = 'green',
            'background data, new basis' = 'orange',
            'background data, old basis' = 'blue',
            'mixed data, old basis' = 'purple',
            'actual background' = 'black',
            'proposed background' = 'red',
            'true model' = 'pink')
plt_T1 <- ggplot(mapping=aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys1[2,]),
                          color = 'mixed data, new basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = bkml_new1,
                          color = 'background data, new basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = bkml_old1,
                          color = 'background data, old basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = bkml_old_full1,
                          color = 'mixed data, old basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                 a = 0, b = 1),
                          color = 'actual background')
            , lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dunif(xs),
                          color = 'proposed background'), lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) + 
  labs(x = "background data",
       y = "",
       color = "Legend") +
  scale_color_manual(values = colors) + 
  ggtitle(bquote(T[1]))
  
#-------------------------------------------------------------------------------

res2 <- constrOptim(theta = rep(0.1,3),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))

back_mle2_new <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,2))


back_mle2_old <- nlm(f = neg_loglikelihood_back, data = data_back, 
                     p = rep(0.01,2))
back_mle2_old_full <- nlm(f = neg_loglikelihood_back, data = data, 
                          p = rep(0.01,2))
(eta <- res2$par[1])
(beta <- res2$par[-1])
(beta_back_new <- back_mle2_new$estimate)
(beta_back_old <- back_mle2_old$estimate)
(beta_back_old_full <- back_mle2_old_full$estimate)

ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml_new2 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back_new, mean = 0.4, sd = 0.1)
})

bkml_old2 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old)
})

bkml_old_full2 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old_full)
})

plt_T2 <- ggplot(mapping=aes(x = data_back)) + 
    geom_histogram(mapping = aes(y = after_stat(density)),
                   bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = xs, y = as.numeric(ys2[2,]),
                          color = 'mixed data, new basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = bkml_new2,
                          color = 'background data, new basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = bkml_old2,
                          color = 'background data, old basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = bkml_old_full2,
                          color = 'mixed data, old basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                 a = 0, b = 1),
                          color = 'actual background')
            , lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dunif(xs),
                          color = 'proposed background'), lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) +
  labs(x = "background data",
       y = "",
       color = "Legend") +
  scale_color_manual(values = colors) +
  ggtitle(bquote(T[1] + T[2]))

#-------------------------------------------------------------------------------

res3 <- constrOptim(theta = rep(0.1,4),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))

back_mle3_new <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,3))

back_mle3_old <- nlm(f = neg_loglikelihood_back, data = data_back, 
                     p = rep(0.01,3))
back_mle3_old_full <- nlm(f = neg_loglikelihood_back, data = data, 
                          p = rep(0.01,3))



(eta <- res3$par[1])
(beta <- res3$par[-1])
(beta_back_new <- back_mle3_new$estimate)
(beta_back_old <- back_mle3_old$estimate)
(beta_back_old_full <- back_mle3_old_full$estimate)

hist(data_back, probability = TRUE, breaks = 20)
xs <- seq(0,1,0.01)

ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml_new3 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back_new, mean = 0.4, sd = 0.1)
})

bkml_old3 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old)
})

bkml_old_full3 <- sapply(xs, function(t)
{
  mod_back(t, beta = beta_back_old_full)
})

plt_T3 <- ggplot(mapping=aes(x = data_back)) + 
    geom_histogram(mapping = aes(y = after_stat(density)),
                   bins = 20, fill = 'steelblue', col = 'black') + 
    geom_line(mapping = aes(x = xs, y = as.numeric(ys3[2,]),
                            color = 'mixed data, new basis'),
              linetype = 'dashed', lwd = 1.1) +
    geom_line(mapping = aes(x = xs, y = bkml_new3,
                            color = 'background data, new basis'),
              linetype = 'dashed', lwd = 1.1) + 
    geom_line(mapping = aes(x = xs, y = bkml_old3,
                            color = 'background data, old basis'),
              linetype = 'dashed', lwd = 1.1) +
    geom_line(mapping = aes(x = xs, y = bkml_old_full3,
                            color = 'mixed data, old basis'),
              linetype = 'dashed', lwd = 1.1) + 
    geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                   a = 0, b = 1),
                            color = 'actual background')
              , lwd = 1.1) + 
    geom_line(mapping = aes(x = xs, y = dunif(xs),
                            color = 'proposed background'), lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) +
    labs(x = "background data",
         y = "",
         color = "Legend") +
    scale_color_manual(values = colors) + ggtitle(bquote(T[1] + T[2] + T[3]))

ggpubr::ggarrange(plt_T1, plt_T2, plt_T3,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = 'bottom')
