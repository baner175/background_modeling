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

mme_estimation <- function(data, nfun, ...){
  lambda <- sapply(data, function(t){
    S(t, ...)$s1
  })|> mean()
  tau <- NULL
  for(j in 1:nfun)
  {
    Tj <- basis[[j+1]]
    coef <- sapply(data, function(t)
      {
      Tj(t, ...)
    })|> mean()
    tau <- c(tau, coef)
  }
  mixture <- function(x, ...)
  {
    val <- 1 + lambda*S(x,...)$s1
    for(i in 1:nfun)
    {
      val <- val+tau[i]*basis[[1+i]](x, ...)
    }
    return(val)
  }
  return(list(coefs = c(lambda, tau), signal = mixture))
}

res1 <- mme_estimation(data, nfun = 1, mean = mean_sig, sd = sd_sig)
res2 <- mme_estimation(data, nfun = 2, mean = mean_sig, sd = sd_sig)
res3 <- mme_estimation(data, nfun = 3, mean = mean_sig, sd = sd_sig)


y1 <- sapply(xs, function(t){
  res1$signal(t, mean = mean_sig, sd = sd_sig)
})
y2 <- sapply(xs, function(t){
  res2$signal(t, mean = mean_sig, sd = sd_sig)
})
y3 <- sapply(xs, function(t){
  res3$signal(t, mean = mean_sig, sd = sd_sig)
})


plt <- ggplot(mapping = aes(x = data)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = xs,
                          y = y_real,
                          color = 'Actual'),
            lwd = 1.1) + 
  geom_line(mapping = aes(x = xs,
                          y = dunif(xs),
                          color = 'proposed bkg'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y1,
                          color = 'T1'),
            lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = y2,
                          color = 'T1 + T2'),
            lwd = 1.1) +
  geom_line(mapping = aes(x = xs, y = y3,
                          color = 'T1 + T2 + T3'),
            lwd = 1.1) +
  labs(x = "background data",
       y = "",
       color = "Legend") + 
  scale_color_manual(values = c('Actual' = 'pink',
                                'proposed bkg' = 'red',
                                'T1' = 'skyblue',
                                'T1 + T2' = 'cyan',
                                'T1 + T2 + T3' = 'darkblue'))
