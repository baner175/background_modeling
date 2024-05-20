rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)

l <- 1; u <- 5

mean_sig <- 3.5; sd_sig <- 0.25
mean_back <- -0.5; sd_back <- 3.25
eta_true <- 0.05

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

f <- function(x)
{
  return(eta_true*dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean_sig, sd = sd_sig) + (1-eta_true)*dtrunc(x, spec = 'norm', a = l, b = u,
                                               mean = mean_back,
                                               sd = sd_back))
}

fb <- function(x)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean_back,
                sd = sd_back))
}

set.seed(12345)
n <- 5e3

back <- rtrunc(n, spec = 'norm',
               mean = mean_back, sd = sd_back,
               a = l, b = u)
sig <- rtrunc(n, spec = 'norm',
              mean = mean_sig, sd = sd_sig,
              a = l, b = u)
u_mask <- rbinom(n, size = 1, prob = eta_true)

obs <- ifelse(u_mask, sig, back)
bw <- 0.2
kde <- kdensity(obs, bw = bw)

fun <- function(t) kde(t) - fs(t)

i1 <- uniroot(fun, interval = c(mean_sig - 3*sd_sig, mean_sig),
              maxiter = 1500)$root
i2 <- uniroot(fun, interval = c(mean_sig, mean_sig + 3*sd_sig),
              maxiter = 1500)$root

curve(kde, l, u, col = 'purple', ylim = c(0,0.4))
curve(fs, l, u, col = 'skyblue', add = TRUE)
abline(lty = 2, v = c(i1, i2))
v1 <- fs(i1)
v2 <- fs(i2)


# lam1 <- log(v1)/(-i1)
# lam2 <- log(v2)/(-i2)
# shift <- 0.1
# a <- weighted.mean(c(i1,i2), w = sqrt(c(v2-shift, v1-shift)))
# c <- (v1-shift)/((i1-a)^2)
# 
# gb <- function(x)
# {
#   val1 <- exp(-lam1*x)
#   val2 <- c*(x-a)^2 + shift
#   val3 <- exp(-lam2*x)
#   inI1 <- l<=x & x < i1
#   inI2 <- i1<=x & x < i2
#   inI3 <- i2<=x & x <= u
#   return(val1*inI1 + val2*inI2 + val3*inI3)
# }
# curve(kde, l, u, col = 'purple', lwd = 2, ylim = c(0,0.6))
# curve(fs, l, u, col = 'skyblue', add = TRUE, lwd = 2)
# curve(gb, l, u, col = 'red', add = TRUE, lwd = 2)
# curve(f, l, u, col = 'black', add = TRUE, lwd = 2)
# curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2)
# abline(lty = 2, v = c(i1, i2))


# shift <- v2
# a <- weighted.mean(c(i1,i2), w = sqrt(c(v2-shift, v1-shift)))
# c <- (v1-shift)/((i1-a)^2)
# mid_quad <- function(x)
# {
#   return(c*(x-a)^2+shift)
# }
# # A2 <- integrate(mid_quad, i1,i2)$value
# A2 <- integrate(mid_quad, i1,u)$value
# # A3 <- v2*(u-i2)
# A3 <- (0.2)*(1-A2)
# # A1 <- 1-A2-A3
# A1 <- 1-A2
# # A1 <- (0.9)*(1-A2)# 0.7293089
# a_lin <- (v1 - A1/(i1-l))*(2/(i1-l))
# b_lin <- v1 - a_lin*i1
# m_lin <- (A3/(u-i2) - v2)*(2/(u-i2))
# z_lin <- v2-m_lin*i2
# gb_lin <- function(x)
# {
#   val1 <- a_lin*x+b_lin
#   val2 <- c*(x-a)^2 + shift
#   val3 <- m_lin*x+z_lin
#   # val3 <- v2
#   #val3 <- exp(-lam2*x)
#   inI1 <- l<=x & x < i1
#   inI2 <- i1<=x & x < i2
#   inI3 <- i2<=x & x <= u
#   return(val1*inI1 + val2*inI2 + val2*inI3)
# }
# 
# integrate(gb_lin, l, u)
# 
# curve(kde, l, u, col = 'purple', lwd = 2, ylim = c(0,0.6), main = 'with gb_lin')
# curve(fs, l, u, col = 'skyblue', add = TRUE, lwd = 2)
# curve(gb_lin, l, u, col = 'red', add = TRUE, lwd = 2)
# curve(f, l, u, col = 'black', add = TRUE, lwd = 2)
# curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2)
# abline(lty = 2, v = c(i1, i2))
# 
# (I1 <- integrate(function(t) (fs(t)/gb_lin(t) - 1)*f(t), l, i1)$value)
# (I2 <- integrate(function(t) (fs(t)/gb_lin(t) - 1)*f(t), i1,i2)$value)
# (I3 <- integrate(function(t) (fs(t)/gb_lin(t) - 1)*f(t), i2,u)$value)
# 
# abs(I1) + abs(I3)
# 
# I2
# 
# #-------------------------------------------------------------------

(I1 <- integrate(function(t) (fs(t)/kde(t) - 1)*f(t), l, i1)$value)
(I2 <- integrate(function(t) (fs(t)/kde(t) - 1)*f(t), i1,i2)$value)
(I3 <- integrate(function(t) (fs(t)/kde(t) - 1)*f(t), i2,u)$value)

abs(I1) + abs(I3)

I2

#-------------------------------------------------------------------
# A <- integrate(kde,l,i2)$value
# A3 <- 1-A
# m_lin <- (A3/(u-i2) - v2)*(2/(u-i2))
# z_lin <- v2-m_lin*i2
# 
# gb_new <- function(x)
# {
#   
#   val1 <- kde(x)
#   val3 <- m_lin*x + z_lin
#   # val3 <- kde(x)
#   inI1UI2 <- l<=x & x < i2
#   inI3 <- i2<=x & x <= u
#   return(val1*inI1UI2 + val3*inI3)
# }
# 
# integrate(gb_new,l,u)
# 
# curve(kde, l, u, col = 'purple', lwd = 2, ylim = c(0,0.6), main = 'with gb_lin')
# curve(fs, l, u, col = 'skyblue', add = TRUE, lwd = 2)
# curve(gb_new, l, u, col = 'red', add = TRUE, lwd = 2)
# curve(f, l, u, col = 'black', add = TRUE, lwd = 2)
# curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2)
# abline(lty = 2, v = c(i1, i2))
# 
# (I1 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, i1)$value)
# (I2 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i1,i2)$value)
# (I3 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i2,u)$value)
# 
# abs(I1) + abs(I3)
# 
# I2
# 
# integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, u)
# 
# S <- function(t) (fs(t)/gb_new(t) - 1)
# 
# mean(sapply(obs, S))

#-------------------------------------------------------------------
p1 <- 1.5
v_p1 <- kde(p1)
A1 <- (p1-l)*v_p1
A2 <- integrate(kde, p1, i2)$value
A <- A1+A2
A3 <- 1-A
m_3 <- (A3/(u-i2) - v2)*(2/(u-i2))
z_3 <- v2-m_3*i2
gb_new <- function(x)
{
  
  val1 <- v_p1
  val2 <- kde(x)
  val3 <- m_3*x + z_3
  inI1 <- x<p1 & x>= l
  inI2 <- x<i2 & x>= p1
  inI3 <- i2<=x & x <= u
  return(val1*inI1 + val2*inI2 + val3*inI3)
}

integrate(gb_new,l,u)

curve(kde, l, u, col = 'purple', lwd = 2, ylim = c(0,0.6), main = 'with gb_lin')
curve(fs, l, u, col = 'skyblue', add = TRUE, lwd = 2)
curve(gb_new, l, u, col = 'red', add = TRUE, lwd = 2)
curve(f, l, u, col = 'black', add = TRUE, lwd = 2)
curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2)
abline(lty = 2, v = c(i1, i2))

(I1 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, i1)$value)
(I2 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i1,i2)$value)
(I3 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i2,u)$value)

abs(I1) + abs(I3)

I2

integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, u)

S <- function(t) (fs(t)/gb_new(t) - 1)

mean(sapply(obs, S))
#-----------------------------------------------------------------
p1 <- 1.5
v_p1 <- kde(p1)
A1 <- (p1-l)*v_p1
A2 <- integrate(kde, p1, i1)$value
shift <- v2
a <- weighted.mean(c(i1,i2), w = sqrt(c(v2-shift, v1-shift)))
c <- (v1-shift)/((i1-a)^2)
mid_quad <- function(x)
{
  return(c*(x-a)^2+shift)
}
sl <- (v2-v1)/(i2-i1)
int <- v1 - sl*i1
rev_quad <- function(x)
{
  diff <- int+sl*x - mid_quad(x)
  return(int+sl*x+diff)
}
A3 <- integrate(rev_quad, i1, i2)$value
A <- A1+A2+A3
A4 <- 1-A
m_3 <- (A4/(u-i2) - v2)*(2/(u-i2))
z_3 <- v2-m_3*i2

gb_new <- function(x)
{
  val1 <- v_p1
  val2 <- kde(x)
  val3 <- rev_quad(x)
  val4 <- m_3*x+z_3
  inI1 <- x<p1 & x>=l
  inI2 <- x<i1 & x>= p1
  inI3 <- x<i2 & x>= i1
  inI4 <- i2<=x & x <= u
  return(val1*inI1 + val2*inI2 + val3*inI3 + val4*inI4)
}

integrate(gb_new,l,u)

curve(kde, l, u, col = 'purple', lwd = 2, ylim = c(0,0.6), main = 'with gb_lin')
curve(fs, l, u, col = 'skyblue', add = TRUE, lwd = 2)
curve(gb_new, l, u, col = 'red', add = TRUE, lwd = 2)
curve(f, l, u, col = 'black', add = TRUE, lwd = 2)
curve(fb, l, u, col = 'brown', add = TRUE, lwd = 2)
abline(lty = 2, v = c(i1, i2))

(I1 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, i1)$value)
(I2 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i1,i2)$value)
(I3 <- integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), i2,u)$value)

abs(I1) + abs(I3)

I2

integrate(function(t) (fs(t)/gb_new(t) - 1)*f(t), l, u)

S <- function(t) (fs(t)/gb_new(t) - 1)
norm_S <- integrate(function(t) (S(t)^2)*gb_new(t),l,u)$value |> sqrt()

(theta <- mean(sapply(obs, S)))

theta/norm_S
