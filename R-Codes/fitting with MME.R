# estimate the mixture and all with MME

source('functions.R')

library(ggplot2)

mean_sig = 0.4
sd_sig = 0.1
data <- make_data(mean_sig = 0.4, 
                  sd_sig = 0.1,
                  mean_back = 5,
                  sd_back = 5,
                  bkg_prop = 0.97,
                  seed = 12345)

data_back <- data$background
data <- data$observed

basis <- c(S,T1,T2,T3)

xs <- seq(0,1,0.01)
y_real <- sapply(xs, actual, mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97)

mme_estimation <- function(data, nfun, mean, sd)
  {
  lambda <- sapply(data, function(t){
    S(t, mean = mean, sd = sd)$s1
  })|> mean()
  norm <- (1/(2*sd*sqrt(pi)))*(pnorm(1,mean, sd/sqrt(2)) - 
                                 pnorm(0,mean, sd/sqrt(2)))/
    (pnorm(1, mean, sd) - pnorm(0, mean, sd))^2 - 1
  norm <- sqrt(norm)
  eta <- lambda/norm
  tau <- NULL
  for(j in 1:nfun)
  {
    Tj <- basis[[j+1]]
    coef <- sapply(data, function(t)
      {
      Tj(t, mean = mean, sd = sd)
    })|> mean()
    tau <- c(tau, coef)
  }
  beta <- tau/(1-eta)
  signal <- function(x)
  {
    mix <- 1 + lambda*S(x, mean = mean, sd = sd)$s1
    bkg <- 1
    for(i in 1:nfun)
    {
      Ti_val <- basis[[1+i]](x, mean = mean, sd = sd)
      mix <- mix+tau[i]*Ti_val
      bkg <- bkg+beta[i]*Ti_val
    }
    
    return(list(mixture = mix, background = bkg))
  }
  return(list(coefs = c(lambda, tau), sig_prop = eta, 
              signal = signal))
}

res1_MME <- mme_estimation(data, nfun = 1, mean = mean_sig, sd = sd_sig)
res1_MLE <- constrOptim(theta = rep(0.1,2),f = neg_loglikelihood,
                        data = data, mean_sig = 0.4, sd_sig = 0.1,
                        grad = score, 
                        ui = rbind(c(1,0),c(-1,0)),
                        ci = c(0,-1))
(eta <- res1_MLE$par[1])
(beta <- res1_MLE$par[-1])

y1_MME <- res1_MME$signal(xs)
ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})
y1_MLE_mix <- as.numeric(ys1[1,])
y1_MLE_back <- as.numeric(ys1[2,])

hs <- hist(data, probability = TRUE, breaks = 20)

plt1 <- ggplot(mapping = aes(x = data)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = y_real,
                          color = 'Actual mixture'),
            lwd = 1.1) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y1_MME$mixture,
                          color = 'MME-mixture'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y1_MLE_mix,
                          color = 'MLE-mixture'),
            lwd = 1.1, linetype = 'longdash') +
  geom_line(mapping = aes(x = xs, y = y1_MME$background,
                          color = 'MME-background'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y1_MLE_back,
                          color = 'MLE-background'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual mixture' = 'black',
                                'Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME-mixture' = 'skyblue',
                                'MLE-mixture' = 'cyan',
                                'MME-background' = 'orange',
                                'MLE-background' = 'pink'))+ 
  ggtitle(bquote(T[1]))
  
#-------------------------------------------------------------------------------
res2_MME <- mme_estimation(data, nfun = 2, mean = mean_sig, sd = sd_sig)
res2_MLE <- constrOptim(theta = rep(0.1,3),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))
(eta <- res2_MLE$par[1])
(beta <- res2_MLE$par[-1])

ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})
y2_MLE_mix <- as.numeric(ys2[1,])
y2_MLE_back <- as.numeric(ys2[2,])

y2_MME <- res2_MME$signal(xs)

plt2 <- ggplot(mapping = aes(x = data)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = y_real,
                          color = 'Actual mixture'),
            lwd = 1.1) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y2_MME$mixture,
                          color = 'MME-mixture'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y2_MLE_mix,
                          color = 'MLE-mixture'),
            lwd = 1.1, linetype = 'longdash') +
  geom_line(mapping = aes(x = xs, y = y2_MME$background,
                          color = 'MME-background'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y2_MLE_back,
                          color = 'MLE-background'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual mixture' = 'black',
                                'Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME-mixture' = 'skyblue',
                                'MLE-mixture' = 'cyan',
                                'MME-background' = 'orange',
                                'MLE-background' = 'pink'))+ 
  ggtitle(bquote(T[1] + T[2]))

#-------------------------------------------------------------------------------
res3_MME <- mme_estimation(data, nfun = 3, mean = mean_sig, sd = sd_sig)

y3_MME <-res3_MME$signal(xs)

res3_MLE <- constrOptim(theta = rep(0.1,4),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))
(eta <- res3_MLE$par[1])
(beta <- res3_MLE$par[-1])
ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})
y3_MLE_mix <- as.numeric(ys3[1,])
y3_MLE_back <- as.numeric(ys3[2,])

plt3 <- ggplot(mapping = aes(x = data)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = y_real,
                          color = 'Actual mixture'),
            lwd = 1.1) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y3_MME$mixture,
                          color = 'MME-mixture'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y3_MLE_mix,
                          color = 'MLE-mixture'),
            lwd = 1.1, linetype = 'longdash') +
  geom_line(mapping = aes(x = xs, y = y3_MME$background,
                          color = 'MME-background'),
            lwd = 1.1, linetype = 'dashed') + 
  geom_line(mapping = aes(x = xs, y = y3_MLE_back,
                          color = 'MLE-background'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual mixture' = 'black',
                                'Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME-mixture' = 'skyblue',
                                'MLE-mixture' = 'cyan',
                                'MME-background' = 'orange',
                                'MLE-background' = 'pink'))+ 
  ggtitle(bquote(T[1] + T[2] + T[3]))


ggpubr::ggarrange(plt1, plt2, plt3,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = 'bottom')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## Analysis on estimating only f_b
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

back_mle_1 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = 0.01)
(beta_back <- back_mle_1$estimate)
bkml_1 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})

bkmm_1 <- mme_estimation(data = data_back, nfun = 1, mean = mean_sig, sd =sd_sig)
backmm_1 <- bkmm_1$signal(xs)$mixture - bkmm_1$coefs[1]*S(xs,
                                                          mean = mean_sig,
                                                          sd = sd_sig)$s1

hs_bkg <- hist(data_back, probability = TRUE, breaks = 20)

plt1_bkg <- ggplot(mapping = aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs_bkg$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = backmm_1,
                          color = 'MME from background only data'),
            lwd = 1.1, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = bkml_1,
                          color = 'MLE from background only data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y1_MME$background,
                          color = 'MME from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y1_MLE_back,
                          color = 'MLE from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME from mixture data' = 'skyblue',
                                'MLE from mixture data' = 'cyan',
                                'MME from background only data' = 'orange',
                                'MLE from background only data' = 'pink'))+ 
  ggtitle(bquote(T[1]))


back_mle_2 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,2))
(beta_back <- back_mle_2$estimate)
bkml_2 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})
bkmm_2 <- mme_estimation(data = data_back, nfun = 2, mean = mean_sig, sd =sd_sig)
backmm_2 <- bkmm_2$signal(xs)$mixture - bkmm_2$coefs[1]*S(xs,
                                                          mean = mean_sig,
                                                          sd = sd_sig)$s1
plt2_bkg <- ggplot(mapping = aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs_bkg$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = backmm_2,
                          color = 'MME from background only data'),
            lwd = 1.1, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = bkml_2,
                          color = 'MLE from background only data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y2_MME$background,
                          color = 'MME from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y2_MLE_back,
                          color = 'MLE from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME from mixture data' = 'skyblue',
                                'MLE from mixture data' = 'cyan',
                                'MME from background only data' = 'orange',
                                'MLE from background only data' = 'pink'))+ 
  ggtitle(bquote(T[1] + T[2]))


back_mle_3 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,3))
(beta_back <- back_mle_3$estimate)
bkml_3 <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})
bkmm_3 <- mme_estimation(data = data_back, nfun = 3, mean = mean_sig, sd = sd_sig)
backmm_3 <- bkmm_3$signal(xs)$mixture - bkmm_3$coefs[1]*S(xs,
                                                          mean = mean_sig,
                                                          sd = sd_sig)$s1
plt3_bkg <- ggplot(mapping = aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black',
                 breaks = hs_bkg$breaks) + 
  geom_line(mapping = aes(x = xs,
                          y = dtruncnorm(xs, a= 0, b= 1,
                                         mean = 5, sd = 5),
                          color = 'Actual background'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = backmm_3,
                          color = 'MME from background only data'),
            lwd = 1.1, linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = bkml_3,
                          color = 'MLE from background only data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y3_MME$background,
                          color = 'MME from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  geom_line(mapping = aes(x = xs, y = y3_MLE_back,
                          color = 'MLE from mixture data'),
            lwd = 1.1, linetype = 'longdash') + 
  labs(x = "mixed data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual background' = 'green',
                                'proposed bkg' = 'red',
                                'MME from mixture data' = 'skyblue',
                                'MLE from mixture data' = 'cyan',
                                'MME from background only data' = 'orange',
                                'MLE from background only data' = 'pink'))+ 
  ggtitle(bquote(T[1] + T[2] + T[3]))

ggpubr::ggarrange(plt1_bkg, plt2_bkg, plt3_bkg,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = 'bottom')
