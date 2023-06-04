
T1 <- function(x, ...)
{
  inner <- integrate(function(t) S(t, ...)$s1*t, 0, 1)$value
  x - inner*S(x, ...)$s1 - 0.5
}

norm1 <- integrate(function(t) T1(t, mean = mean_sig,
                                  sd = sd_sig)^2, 0, 1)$value
norm1 <- sqrt(norm1)


T2 <- function(x, ...)
{
  inner1 <- integrate(function(t) (T1(t, ...)/norm1)*(t^2), 0,1)$value
  inner2 <- integrate(function(t) S(t, ...)$s1*(t^2), 0,1)$value
  
  x^2 - T1(x,...)*inner1/norm1 - S(x,...)$s1*inner2 - (1/3)
}

norm2 <- integrate(function(t) T2(t, mean = mean_sig,
                                  sd = sd_sig)^2, 0, 1)$value
norm2 <- sqrt(norm2)

T3 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^3)/norm1, 0,1)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^3)/norm2, 0,1)$value
  inner3 <- integrate(function(t) S(t, ...)$s1*(t^3), 0,1)$value
  
  x^3 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - S(x,...)$s1*inner3 - (1/4)
}

norm3 <- integrate(function(t) T3(t, mean = mean_sig,
                                  sd = sd_sig)^2, 0, 1)$value
norm3 <- sqrt(norm3)


mod <- function(x, eta, beta, ...)
{
  blen <- length(beta)
  if(blen == 1) Tvec <- T1(x,...)/norm1
  if(blen == 2) Tvec <- c(T1(x,...)/norm1, T2(x,...)/norm2)
  if(blen == 3) Tvec <- c(T1(x,...)/norm1, T2(x,...)/norm2,T3(x,...)/norm3)
  
  mix <- 1 + eta*S(x, ...)$s + (1-eta)*as.numeric(crossprod(beta,Tvec))
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

  blen <- length(beta)
  if(blen == 1){
    Tmat <- sapply(data, function(t) T1(t, mean = mean_sig, sd = sd_sig)/norm1)
    dim(Tmat) <- c(1, length(data))
    der_beta <- (1-eta)*Tmat%*%rec

  }
  if(blen == 2){
    Tmat <- rbind(sapply(data, function(t) T1(t, mean = mean_sig, sd = sd_sig)/norm1),
                  sapply(data, function(t) T2(t, mean = mean_sig, sd = sd_sig)/norm2))
    der_beta <- (1-eta)*Tmat%*%rec

  }
  if(blen == 3){
    Tmat <- rbind(sapply(data, function(t) T1(t, mean = mean_sig, sd = sd_sig)/norm1),
                  sapply(data, function(t) T2(t, mean = mean_sig, sd = sd_sig)/norm2),
                  sapply(data, function(t) T3(t, mean = mean_sig, sd = sd_sig)/norm3))
    
    der_beta <- (1-eta)*Tmat%*%rec
  }

  der_eta_num <- sapply(data, function(t)
    {
      index <- which(data == t)
      Tvec <- Tmat[,index]
      return(S(t, mean = mean_sig, sd = sd_sig)$s -
      as.numeric(crossprod(beta, Tvec)))
  })

  return(-1*c(
    sum(der_eta_num*rec), as.vector(der_beta)
  ))
}


mod_back_new <- function(x, beta, ...)
{
  blen <- length(beta)
  if(blen == 1) Tvec <- T1(x,...)/norm1
  if(blen == 2) Tvec <- c(T1(x,...)/norm1, T2(x,...)/norm2)
  if(blen == 3) Tvec <- c(T1(x,...)/norm1, T2(x,...)/norm2,T3(x,...)/norm3)
  
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

