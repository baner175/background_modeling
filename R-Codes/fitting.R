# source('functions.R')

par(mfrow = c(1,3))


data <- make_data(mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97,
                 seed = 12345)

data_back <- data$background
data <- data$observed


# hist(data, probability = TRUE, breaks = 20)

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

back_mle_1 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = 0.01)

# old_basis_mix_1 <- nlm(f = neg_loglikelihood_back,
#                      data = data, p = 0.01)
  
  # optim(par = 0.1, 
  #       fn = neg_loglikelihood_back,
  #       data =  data_back,
  #       method = 'Nelder-Mead')
  
(eta <- res1$par[1])
(beta <- res1$par[-1])
(beta_back <- back_mle_1$estimate)
# (beta_back_mix <- old_basis_mix_1$estimate)

hist(data, probability = TRUE, breaks = 20)

ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml <- sapply(xs, function(t)
  {
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})

# mix_old <- sapply(xs, function(t)
#   {
#   mod_back(t, beta = beta_back_mix)
# })
lines(xs, ys1[1,], col = 'blue', lwd =2) # estimated mixture
lines(xs, ys1[2,], col = 'green', lwd =2, lty = 2) # estimated background from the mixed data
lines(xs, bkml, col = 'orange', lwd =2, lty = 2) # estimated background from background data
lines(xs, y_real, col = 'black', lwd = 2, lty = 1) # plotting the actual mixture signal
lines(xs, dtruncnorm(xs, mean = 5, sd = 5,
                     a = 0, b = 1),
      col = 'purple', lwd = 2, lty = 1) # plotting the actual background
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1) # proposed background
#-------------------------------------------------------------------------------

res2 <- constrOptim(theta = rep(0.1,3),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))

back_mle_2 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,2))


(eta <- res2$par[1])
(beta <- res2$par[-1])
(beta_back <- back_mle_2$estimate)

hist(data, probability = TRUE, breaks = 20)
xs <- seq(0,1,0.01)

ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})

lines(xs, ys2[1,], col = 'blue', lwd =2)
lines(xs, ys2[2,], col = 'green', lwd =2, lty = 2)
lines(xs, bkml, col = 'orange', lwd =2, lty = 2)
lines(xs, y_real, col = 'black', lwd = 2, lty = 1)
lines(xs, dtruncnorm(xs, mean = 5, sd = 5,
                     a = 0, b = 1),
      col = 'purple', lwd = 2, lty = 1)
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1)
#-------------------------------------------------------------------------------

res3 <- constrOptim(theta = rep(0.1,4),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))

back_mle_3 <- nlm(f = neg_loglikelihood_back_new, data = data_back,
                  mean = 0.4, sd_sig = 0.1, p = rep(0.01,3))


(eta <- res3$par[1])
(beta <- res3$par[-1])
(beta_back <- back_mle_3$estimate)
hist(data, probability = TRUE, breaks = 20)
xs <- seq(0,1,0.01)

ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

bkml <- sapply(xs, function(t)
{
  mod_back_new(t, beta = beta_back, mean = 0.4, sd = 0.1)
})

lines(xs, ys3[1,], col = 'blue', lwd =2)
lines(xs, ys3[2,], col = 'green', lwd =2, lty = 2)
lines(xs, bkml, col = 'orange', lwd =2, lty = 2)
lines(xs, y_real, col = 'black', lwd = 2, lty = 1)
lines(xs, dtruncnorm(xs, mean = 5, sd = 5,
                     a = 0, b = 1),
      col = 'purple', lwd = 2, lty = 1)
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1)
# lines(xs, sapply(xs, mod, eta = eta, beta = beta, mean = 0.4, sd = 0.1),
#       col = 'orange', lwd = 2)
