library(truncnorm)

mean_back <- 5; sd_back <- 5
mean_sig <- 0.4; sd_sig <- 0.1
eta_true <- 0.03

l <- 0; u <- 1
#-------------------------------------------------------------------------------
# computing the norm numerically; can be replaced with theoretical value later
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
  inner <- integrate(function(t) S1(t, ...)*t, l, u)$value
  inner_last <- integrate(function(t) t, l,u)$value
  x - inner*S1(x, ...) - inner_last
}
norm1 <- integrate(function(t) T1(t, mean = mean_sig,
                                  sd = sd_sig)^2, l, u)$value
norm1 <- sqrt(norm1)

T2 <- function(x, ...)
{
  inner1 <- integrate(function(t) (T1(t, ...)/norm1)*(t^2), l,u)$value
  inner2 <- integrate(function(t) S1(t, ...)*(t^2), l,u)$value
  inner_last <- integrate(function(t) t^2, l,u)$value
  
  x^2 - inner1*T1(x,...)/norm1 - inner2*S1(x, ...)  - inner_last
}

norm2 <- integrate(function(t) T2(t, mean = mean_sig,
                                  sd = sd_sig)^2, l, u)$value
norm2 <- sqrt(norm2)

T3 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^3)/norm1, l,u)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^3)/norm2, l,u)$value
  inner3 <- integrate(function(t) S1(t, ...)*(t^3), l,u)$value
  inner_last <- integrate(function(t) t^3, l,u)$value
  
  x^3 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - S1(x,...)*inner3 - 
    inner_last
}

norm3 <- integrate(function(t) T3(t, mean = mean_sig,
                                  sd = sd_sig)^2, l, u)$value
norm3 <- sqrt(norm3)

T4 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^4)/norm1, l,u)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^4)/norm2, l,u)$value
  inner3 <- integrate(function(t) T3(t, ...)*(t^4)/norm3, l,u)$value
  inner4 <- integrate(function(t) S1(t, ...)*(t^4), l,u)$value
  inner_last <- integrate(function(t) t^4, l,u)$value
  
  x^4 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - T3(x,...)*inner3/norm3 -
    S1(x,...)*inner4 - inner_last
}

norm4 <- integrate(function(t) T4(t, mean = mean_sig,
                                  sd = sd_sig)^2, l, u)$value
norm4 <- sqrt(norm4)

T_basis <- c(T1, T2, T3, T4)
norm_vec <- c(norm1, norm2, norm3, norm4)

mod <- function(x, eta, beta, ...)
{
  blen <- length(beta)
  Tvec <- c()
  for(j in 1:blen) {
    Tvec <- c(Tvec,
              T_basis[[j]](x, ...)/norm_vec[j])
  }
  mix <- 1 + eta*norm_S*S1(x, ...)+ (1-eta)*as.numeric(crossprod(beta,Tvec))
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(list(mix = mix, back = back))
}


neg_loglikelihood <- function(theta, data, mean_sig, sd_sig)
{
  eta = theta[1]
  beta = theta[-1]
  y <- sapply(data, function(t) mod(x = t, eta = eta,
                                    beta = beta, mean = mean_sig, 
                                    sd = sd_sig)$mix)
  -sum(log(y))
}

score <- function(theta, data, mean_sig, sd_sig)
{
  eta = theta[1]
  beta = theta[-1]

  rec <- sapply(data, function(t)
  {
    1/mod(x = t,  eta = eta,
          beta = beta, mean = mean_sig,
          sd = sd_sig)$mix
  })

  Tmat <- c()
  blen <- length(beta)
  for(j in 1:blen)
  {
    Tmat <- rbind(Tmat,
                  sapply(data, function(t) T_basis[[j]](t, mean = mean_sig,
                                                        sd = sd_sig)/norm_vec[j]))
  }
  der_beta <- (1-eta)*Tmat%*%rec
  beta_T_dot <- apply(Tmat, 2,
                      function(t) as.numeric(crossprod(beta, t)))
  der_eta_num <- norm_S*S1(data, mean = mean_sig, sd = sd_sig) - beta_T_dot

  return(-1*c(
    sum(der_eta_num*rec), as.vector(der_beta)
  ))
}

mod_back_new <- function(x, beta, ...)
{
  blen <- length(beta)
  Tvec <- c()
  for(j in 1:blen) {
    Tvec <- c(Tvec,
              T_basis[[j]](x, ...)/norm_vec[j])
  }
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(back)
}

neg_loglikelihood_back_new <- function(theta, data, mean_sig, sd_sig)
{
  beta <- theta
  y <- sapply(data, function(t) mod_back_new(x = t,
                                             beta = beta, mean = mean_sig,
                                             sd = sd_sig))
  -sum(log(y))
}


