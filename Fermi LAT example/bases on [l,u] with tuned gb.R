source('designing gb with 2 bumps.R')
# Area under gb
integrate(gb,l,u)

# Function to calculate norm under Gb
calc_norm_gb <- function(fun, ...)
{
  integrate(function(t) gb(t)*fun(t,...)^2, l, u)$value |> sqrt()
}

# Convarting fs into S1:
norm_S <- calc_norm_gb(function(t) fs(t, mean= mean_sig, sd = sd_sig)/gb(t) - 1)

S1 <- function(x, ...)
{
  f_sig <- fs(x, ...)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}

# Basis functions:
T1_inner_S1 <- integrate(function(t) t*S1(t,
                                          mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T1_inner_1 <- integrate(function(t) t*gb(t), l, u)$value
T1 <- function(x)
{
  return(x - T1_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) - T1_inner_1)
}

T1_norm <- calc_norm_gb(T1)
T1_normed <- function(x) T1(x)/T1_norm

T2_inner_S1 <- integrate(function(t) (t^2)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T2_inner_T1 <- integrate(function(t) (t^2)*T1_normed(t)*gb(t), l, u)$value
T2_inner_1 <- integrate(function(t) (t^2)*gb(t), l, u)$value

T2 <- function(x)
{
  return(x^2 - T2_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T2_inner_T1*T1_normed(x) - T2_inner_1)
}

T2_norm <- calc_norm_gb(T2)
T2_normed <- function(x) T2(x)/T2_norm


T3_inner_S1 <- integrate(function(t) (t^3)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T3_inner_T1 <- integrate(function(t) (t^3)*T1_normed(t)*gb(t), l, u)$value
T3_inner_T2 <- integrate(function(t) (t^3)*T2_normed(t)*gb(t), l, u)$value
T3_inner_1 <- integrate(function(t) (t^3)*gb(t), l, u)$value

T3 <- function(x)
{
  return(x^3 - T3_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T3_inner_T1*T1_normed(x) - T3_inner_T2*T2_normed(x) -
           T3_inner_1)
}

T3_norm <- calc_norm_gb(T3)
T3_normed <- function(x) T3(x)/T3_norm


T4_inner_S1 <- integrate(function(t) (t^4)*S1(t,
                                              mean = mean_sig, sd = sd_sig)*gb(t), l, u)$value
T4_inner_T1 <- integrate(function(t) (t^4)*T1_normed(t)*gb(t), l, u)$value
T4_inner_T2 <- integrate(function(t) (t^4)*T2_normed(t)*gb(t), l, u)$value
T4_inner_T3 <- integrate(function(t) (t^4)*T3_normed(t)*gb(t), l, u)$value
T4_inner_1 <- integrate(function(t) (t^4)*gb(t), l, u)$value

T4 <- function(x)
{
  return(x^4 - T4_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T4_inner_T1*T1_normed(x) - T4_inner_T2*T2_normed(x) -
           T4_inner_T3*T3_normed(x) - T4_inner_1)
}

T4_norm <- calc_norm_gb(T4)
T4_normed <- function(x) T4(x)/T4_norm


T5_inner_S1 <- integrate(function(t) (t^5)*S1(t,
                                              mean = mean_sig,
                                              sd = sd_sig)*gb(t), l, u)$value
T5_inner_T1 <- integrate(function(t) (t^5)*T1_normed(t)*gb(t), l, u)$value
T5_inner_T2 <- integrate(function(t) (t^5)*T2_normed(t)*gb(t), l, u)$value
T5_inner_T3 <- integrate(function(t) (t^5)*T3_normed(t)*gb(t), l, u)$value
T5_inner_T4 <- integrate(function(t) (t^5)*T4_normed(t)*gb(t), l, u)$value
T5_inner_1 <- integrate(function(t) (t^5)*gb(t), l, u)$value

T5 <- function(x)
{
  return(x^5 - T5_inner_S1*S1(x, mean = mean_sig, sd = sd_sig) -
           T5_inner_T1*T1_normed(x) - T5_inner_T2*T2_normed(x) -
           T5_inner_T3*T3_normed(x) - T5_inner_T4*T4_normed(x) -
           T5_inner_1)
}

T5_norm <- calc_norm_gb(T5)
T5_normed <- function(x) T5(x)/T5_norm

# The basis:
T_basis <- c(T1_normed, T2_normed, T3_normed, T4_normed, T5_normed)

# model:
mod <- function(x, theta, tau)
{
  tau_len <- length(tau)
  Tvec <- c()
  for(j in 1:tau_len) {
    Tvec <- c(Tvec,
              T_basis[[j]](x))
  }
  mix <- 1 + theta*S1(x, mean = mean_sig, sd = sd_sig) + as.numeric(crossprod(tau,Tvec))
  mix <- mix*gb(x)
  eta <- theta/norm_S
  beta <- tau/(1-eta)
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  back <- back*gb(x)
  return(list(mix = mix, back = back))
}

#log-likelihood:
log_lik <- function(data, theta, tau)
{
  y <- sapply(data, function(t) mod(t, theta = theta, tau = tau)$mix)
  return(sum(log(y)))
}

