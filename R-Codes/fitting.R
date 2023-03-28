# source('functions.R')

par(mfrow = c(1,3))


data <- make_data(mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97,
                 seed = 12345)

data <- data$observed


# hist(data, probability = TRUE, breaks = 20)

xs <- seq(0,1,0.01)
y_real <- sapply(xs, actual, mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97)


res <- constrOptim(theta = rep(0.1,2),f = neg_loglikelihood,
                   data = data, mean_sig = 0.4, sd_sig = 0.1,
                   grad = score, 
                   ui = rbind(c(1,0),c(-1,0)),
                   ci = c(0,-1))
(eta <- res$par[1])
(beta <- res$par[-1])

hist(data, probability = TRUE, breaks = 20)

ys1 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

lines(xs, ys1[1,], col = 'blue', lwd =2)
lines(xs, ys1[2,], col = 'green', lwd =2, lty = 2)
lines(xs, y_real, col = 'black', lwd = 2, lty = 1)
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1)
#-------------------------------------------------------------------------------

res2 <- constrOptim(theta = rep(0.1,3),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0),c(-1,0,0)),
                    ci = c(0,-1))

(eta <- res2$par[1])
(beta <- res2$par[-1])

hist(data, probability = TRUE, breaks = 20)
xs <- seq(0,1,0.01)

ys2 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

lines(xs, ys2[1,], col = 'blue', lwd =2)
lines(xs, ys2[2,], col = 'green', lwd =2, lty = 2)
lines(xs, y_real, col = 'black', lwd = 2, lty = 1)
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1)
#-------------------------------------------------------------------------------

res3 <- constrOptim(theta = rep(0.1,4),f = neg_loglikelihood,
                    data = data, mean_sig = 0.4, sd_sig = 0.1,
                    grad = score, 
                    ui = rbind(c(1,0,0,0),c(-1,0,0,0)),
                    ci = c(0,-1))

(eta <- res3$par[1])
(beta <- res3$par[-1])

hist(data, probability = TRUE, breaks = 20)
xs <- seq(0,1,0.01)

ys3 <- sapply(xs, function(t)
{
  mod(t, eta = eta, beta = beta, mean = 0.4, sd = 0.1)
})

lines(xs, ys3[1,], col = 'blue', lwd =2)
lines(xs, ys3[2,], col = 'green', lwd =2, lty = 2)
lines(xs, y_real, col = 'black', lwd = 2, lty = 1)
lines(xs, dunif(xs), col = 'red', lwd = 2, lty = 1)
# lines(xs, sapply(xs, mod, eta = eta, beta = beta, mean = 0.4, sd = 0.1),
#       col = 'orange', lwd = 2)
