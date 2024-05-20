rm(list = ls())
library(truncdist)
library(VGAM)

l <- 1; u <- 5

mean_sig <- 4; sd_sig <- 0.25
mean_back <- -0.5; sd_back <- 3.25
rate_gb <- 0.25
eta_true <- 0.05
#--------------------------------------------------------

gb <- function(x)
{
  return(dtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u))
}

Gb <- function(x)
{
  return(ptrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u))
}

Gb_I <- function(x)
{
  return(qtrunc(x, spec = 'exp', rate = rate_gb, a = l, b = u))
}

fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean, sd = sd))
}

calc_norm <- function(fun, ...)
{
  integrate(function(t) fun(t,...)^2, 0, 1)$value |> sqrt()
}

calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

norm_S <- calc_norm_gb(function(t) fs(t, mean= mean_sig, sd = sd_sig)/gb(t) - 1)

S1 <- function(x, ...)
{
  f_sig <- fs(x, ...)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

T1_inner_S1 <- integrate(function(t) t*S1(Gb_I(t),
                                          mean = mean_sig, sd = sd_sig), 0, 1)$value
T1 <- function(x)
{
  Gb_inv <- Gb_I(x)
  return(x - T1_inner_S1*S1(Gb_inv, mean = mean_sig, sd = sd_sig) - 0.5)
}

T1_norm <- calc_norm(T1)
T1_normed <- function(x) T1(x)/T1_norm

T2_inner_S1 <- integrate(function(t) (t^2)*S1(Gb_I(t),
                                              mean = mean_sig, sd = sd_sig), 0, 1)$value
T2_inner_T1 <- integrate(function(t) (t^2)*T1_normed(t), 0, 1)$value

T2 <- function(x)
{
  Gb_inv <- Gb_I(x)
  return(x^2 - T2_inner_S1*S1(Gb_inv, mean = mean_sig, sd = sd_sig) -
           T2_inner_T1*T1_normed(x) - (1/3))
}

T2_norm <- calc_norm(T2)
T2_normed <- function(x) T2(x)/T2_norm


T3_inner_S1 <- integrate(function(t) (t^3)*S1(Gb_I(t),
                                              mean = mean_sig, sd = sd_sig), 0, 1)$value
T3_inner_T1 <- integrate(function(t) (t^3)*T1_normed(t), 0, 1)$value
T3_inner_T2 <- integrate(function(t) (t^3)*T2_normed(t), 0, 1)$value

T3 <- function(x)
{
  Gb_inv <- Gb_I(x)
  return(x^3 - T3_inner_S1*S1(Gb_inv, mean = mean_sig, sd = sd_sig) -
           T3_inner_T1*T1_normed(x) - T3_inner_T2*T2_normed(x) - 
           (1/4))
}

T3_norm <- calc_norm(T3)
T3_normed <- function(x) T3(x)/T3_norm


T4_inner_S1 <- integrate(function(t) (t^4)*S1(Gb_I(t),
                                              mean = mean_sig, sd = sd_sig), 0, 1)$value
T4_inner_T1 <- integrate(function(t) (t^4)*T1_normed(t), 0, 1)$value
T4_inner_T2 <- integrate(function(t) (t^4)*T2_normed(t), 0, 1)$value
T4_inner_T3 <- integrate(function(t) (t^4)*T3_normed(t), 0, 1)$value

T4 <- function(x)
{
  Gb_inv <- Gb_I(x)
  return(x^4 - T4_inner_S1*S1(Gb_inv, mean = mean_sig, sd = sd_sig) -
           T4_inner_T1*T1_normed(x) - T4_inner_T2*T2_normed(x) - 
           T4_inner_T3*T3_normed(x) - (1/5))
}

T4_norm <- calc_norm(T4)
T4_normed <- function(x) T4(x)/T4_norm


T5_inner_S1 <- integrate(function(t) (t^5)*S1(Gb_I(t),
                                              mean = mean_sig,
                                              sd = sd_sig), 0, 1)$value
T5_inner_T1 <- integrate(function(t) (t^5)*T1_normed(t), 0, 1)$value
T5_inner_T2 <- integrate(function(t) (t^5)*T2_normed(t), 0, 1)$value
T5_inner_T3 <- integrate(function(t) (t^5)*T3_normed(t), 0, 1)$value
T5_inner_T4 <- integrate(function(t) (t^5)*T4_normed(t), 0, 1)$value

T5 <- function(x)
{
  Gb_inv <- Gb_I(x)
  return(x^5 - T5_inner_S1*S1(Gb_inv, mean = mean_sig, sd = sd_sig) -
           T5_inner_T1*T1_normed(x) - T5_inner_T2*T2_normed(x) - 
           T5_inner_T3*T3_normed(x) - T5_inner_T4*T4_normed(x) - (1/6))
}

T5_norm <- calc_norm(T5)
T5_normed <- function(x) T5(x)/T5_norm


T_basis <- c(T1_normed, T2_normed, T3_normed, T4_normed, T5_normed)

mod <- function(x, theta, tau,  ...)
{
  tau_len <- length(tau)
  Tvec <- c()
  for(j in 1:tau_len) {
    Tvec <- c(Tvec,
              T_basis[[j]](Gb(x)))
  }
  mix <- 1 + theta*S1(x, ...) + as.numeric(crossprod(tau,Tvec))
  mix <- mix*gb(x)
  eta <- theta/norm_S
  beta <- tau/(1-eta)
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  back <- back*gb(x)
  return(list(mix = mix, back = back))
}

neg_loglikelihood <- function(params, data = obs, mean_sig = mean_sig,
                              sd_sig = sd_sig)
{
  theta = params[1]
  tau = params[-1]
  y <- sapply(data, function(t) mod(x = t, theta = theta,
                                    tau = tau, 
                                    mean = mean_sig, 
                                    sd = sd_sig)$mix)
  -sum(log(y))
}

  