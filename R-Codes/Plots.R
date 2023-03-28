# source('functions.R')

x <- seq(0,1,0.01)

ys1 <- sapply(x, function(t, mean, sd){
  S(t, mean, sd)$s1
}, mean = 0.4, sd = 0.1)


t1 <- sapply(x, T1, mean = 0.4, sd = 0.1)/norm1
t2 <- sapply(x, T2, mean = 0.4, sd = 0.1)/norm2
t3 <- sapply(x, T3, mean = 0.4, sd = 0.1)/norm3

par(mfrow = c(2,2))
plot(x,ys1,type = 'l', lwd = 2, 
     xlab = bquote('S[1]'))
plot(x,t1,type = 'l', lwd = 2)
plot(x,t2,type = 'l', lwd = 2)
plot(x,t3,type = 'l', lwd = 2)
