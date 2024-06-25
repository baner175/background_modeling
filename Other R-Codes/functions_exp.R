rate_sig <- 2
sig_prop <- 0.05

#-------------------------------------------------------------------------------

library(truncdist)

S <- function(x, rate)
{
  f_sig <- dtrunc(x, spec = 'exp', rate = rate, a = 2, b = 3)
  norm <- (rate/2)*(exp(-2*rate) + 
                      exp(-3*rate))/(exp(-2*rate) - 
                                       exp(-3*rate)) -1
  norm <- sqrt(norm)
  return(list('s' = f_sig - 1, 's1' = (f_sig-1)/norm))
}

T1 <- function(x, ...)
{
  inner <- integrate(function(t) S(t, ...)$s1*t, 2, 3)$value
  x - inner*S(x, ...)$s1 - (5/2)
}

norm1 <- integrate(function(t) T1(t, rate = rate_sig)^2, 2, 3)$value
norm1 <- sqrt(norm1)


T2 <- function(x, ...)
{
  inner1 <- integrate(function(t) (T1(t, ...)/norm1)*(t^2), 2,3)$value
  inner2 <- integrate(function(t) S(t, ...)$s1*(t^2), 2,3)$value
  
  x^2 - T1(x,...)*inner1/norm1 - S(x,...)$s1*inner2 - (19/3)
}

norm2 <- integrate(function(t) T2(t, rate = rate_sig)^2, 2, 3)$value
norm2 <- sqrt(norm2)

T3 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^3)/norm1, 2,3)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^3)/norm2, 2,3)$value
  inner3 <- integrate(function(t) S(t, ...)$s1*(t^3), 2,3)$value
  
  x^3 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - S(x,...)$s1*inner3 - (65/4)
}

norm3 <- integrate(function(t) T3(t,rate = rate_sig)^2, 0, 1)$value
norm3 <- sqrt(norm3)

T4 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^4)/norm1, 2,3)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^4)/norm2, 2,3)$value
  inner3 <- integrate(function(t) T3(t, ...)*(t^4)/norm3, 2,3)$value
  inner4 <- integrate(function(t) S(t, ...)$s1*(t^4), 2,3)$value
  
  x^4 - 
    T1(x,...)*inner1/norm1 - 
    T2(x,...)*inner2/norm2 - 
    T3(x,...)*inner3/norm3 - 
    S(x,...)$s1*inner4 - 211/5
}

norm4 <- integrate(function(t) T4(t, rate = rate_sig)^2, 2, 3)$value
norm4 <- sqrt(norm4)

T_basis <- c(T1, T2, T3, T4)
norm_vec <- c(norm1, norm2, norm3, norm4)


mod <- function(x, eta, beta,  ...)
{
  blen <- length(beta)
  Tvec <- c()
  for(j in 1:blen) {
    Tvec <- c(Tvec,
              T_basis[[j]](x, ...)/norm_vec[j])
  }
  mix <- 1 + eta*S(x, ...)$s + (1-eta)*as.numeric(crossprod(beta,Tvec))
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(list(mix = mix, back = back))
  
}

neg_loglikelihood <- function(theta, data, rate_sig)
{
  eta = theta[1]
  beta = theta[-1]
  y <- sapply(data, function(t) mod(x = t, eta = eta,
                                    beta = beta, 
                                    rate = rate_sig)$mix)
  -sum(log(y))
}

score <- function(theta, data, rate_sig)
{
  eta = theta[1]
  beta = theta[-1]
  
  rec <- sapply(data, function(t)
  {
    1/mod(x = t,  eta = eta,
          beta = beta, 
          rate = rate_sig)$mix
  })
  
  blen <- length(beta)
  Tmat <- c()
  for(j in 1:blen)
  {
    Tmat <- rbind(Tmat,
                  sapply(data, function(t) T_basis[[j]](t, rate = rate_sig)/norm_vec[j]))
  }
  der_beta <- (1-eta)*Tmat%*%rec
  
  beta_T_dot <- apply(Tmat, 2,
                      function(t) as.numeric(crossprod(beta, t)))
  der_eta_num <- S(data, rate = rate_sig)$s - beta_T_dot
  
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

neg_loglikelihood_back_new <- function(theta, data, rate_sig)
{
  beta <- theta
  y <- sapply(data, function(t) mod_back_new(x = t,
                                             beta = beta, 
                                             rate = rate_sig))
  -sum(log(y))
}

