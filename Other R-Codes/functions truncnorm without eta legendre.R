library(truncnorm)

mean_back <- 5; sd_back <- 5
mean_sig <- 0.4; sd_sig <- 0.1
eta_true <- 0.03

l <- 0; u <- 1

P1 <- function(x)
{
  x - 0.5*(u^2 - l^2)
}
P1(1)
(norm_P1 <- integrate(function(t) P1(t)^2, l, u)$value |> sqrt())

P2 <- function(x)
{
  inner_P1 <- integrate(function(t) (t^2)*P1(t)/norm_P1, l, u)$value
  inner_1 <- integrate(function(t) t^2, l, u)$value
  x^2 - P1(x)*inner_P1/norm_P1 - inner_1
}

(norm_P2 <- integrate(function(t) P2(t)^2, l, u)$value |> sqrt())

P3 <- function(x)
{
  inner_1 <- integrate(function(t) t^3, l, u)$value
  inner_P1 <- integrate(function(t) (t^3)*P1(t)/norm_P1, l, u)$value
  inner_P2 <- integrate(function(t) (t^3)*P2(t)/norm_P2, l, u)$value
  
  x^3 - inner_P1*P1(x)/norm_P1 - inner_P2*P2(x)/norm_P2 - inner_1
}
(norm_P3 <- integrate(function(t) P3(t)^2, l, u)$value |> sqrt())

P1_normed <- function(x) P1(x)/norm_P1
P2_normed <- function(x) P2(x)/norm_P2
P3_normed <- function(x) P3(x)/norm_P3

# curve(P1_normed, l, u)
# curve(P2_normed, l, u, add = TRUE, col = 'red')
# curve(P3_normed, l, u, add = TRUE, col = 'blue')

norm_S <- integrate(function(t) {(dtruncnorm(t, mean = mean_sig, 
                                             sd = sd_sig, a = l, b = u) - 1)^2},
                    lower = l, upper = u)$value |> sqrt()


S1 <- function(x, mean, sd)
{
  f_sig <- dtruncnorm(x, mean = mean, sd = sd, a = l, b = u)
  return((f_sig-1)/norm_S)
}

T1 <- function(x, ...)
{
  inner_S1 <- integrate(function(t) S1(t, ...)*P1_normed(t), l, u)$value
  P1_normed(x) - inner_S1*S1(x, ...)
}

norm_T1 <- integrate(function(t) T1(t, mean = mean_sig,
                                    sd = sd_sig)^2, l, u)$value |> sqrt()
T1_normed <- function(x,...)
{
  T1(x,...)/norm_T1
}

T2 <- function(x, ...)
{
  inner_S1 <- integrate(function(t) S1(t, ...)*P2_normed(t), l, u)$value
  inner_T1 <- integrate(function(t) T1_normed(t, ...)*P2_normed(t), l, u)$value
  
  P2_normed(x) - T1_normed(x,...)*inner_T1 - S1(x,...)*inner_S1
}

norm_T2 <- integrate(function(t) T2(t, mean = mean_sig,
                                    sd = sd_sig)^2, l, u)$value |> sqrt()
T2_normed <- function(x,...)
{
  T2(x,...)/norm_T2
}

T3 <- function(x, ...)
{
  inner_S1 <- integrate(function(t) S1(t, ...)*P3_normed(t), l, u)$value
  inner_T1 <- integrate(function(t) T1_normed(t, ...)*P3_normed(t), l, u)$value
  inner_T2 <- integrate(function(t) T2_normed(t, ...)*P3_normed(t), l, u)$value
  
  P3_normed(x) - inner_T1*T1_normed(x,...) - inner_T2*T2_normed(x,...) - inner_S1*S1(x,...)
}

norm_T3 <- integrate(function(t) T3(t, mean = mean_sig,
                                    sd = sd_sig)^2, l, u)$value |> sqrt()
T3_normed <- function(x, ...)
{
  T3(x,...)/norm_T3
}

T_basis <- c(T1_normed, T2_normed, T3_normed)

mod <- function(x, theta, tau,  ...)
{
  tau_len <- length(tau)
  Tvec <- c()
  for(j in 1:tau_len) {
    Tvec <- c(Tvec,
              T_basis[[j]](x, ...))
  }
  mix <- 1 + theta*S1(x, ...) + as.numeric(crossprod(tau,Tvec))
  eta <- theta/norm_S
  beta <- tau/(1-eta)
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(list(mix = mix, back = back))
}

neg_loglikelihood <- function(params, data, mean_sig, sd_sig)
{
  theta = params[1]
  tau = params[-1]
  y <- sapply(data, function(t) mod(x = t, theta = theta,
                                    tau = tau, 
                                    mean = mean_sig,
                                    sd = sd_sig)$mix)
  -sum(log(y))
}

score <- function(params, data, mean_sig, sd_sig)
{
  theta <-  params[1]
  tau <-  params[-1]
  
  rec <- sapply(data, function(t)
  {
    1/mod(x = t,  theta = theta,
          tau = tau, mean = mean_sig,
          sd = sd_sig)$mix
  })
  
  tau_len <- length(tau)
  
  Tmat <- c()
  for(j in 1:tau_len)
  {
    Tmat <- rbind(Tmat,
                  sapply(data, function(t) T_basis[[j]](t, mean = mean_sig,
                                                        sd = sd_sig)))
  }
  der_theta_num <- sapply(data, function(t) S1(t, mean = mean_sig, sd = sd_sig))
  der_beta <- Tmat%*%rec
  return(-1*c(
    sum(der_theta_num*rec), as.vector(der_beta)
  ))
}



