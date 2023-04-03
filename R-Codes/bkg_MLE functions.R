
T1_s <- function(x)
{
  x - 0.5
}
norm1_back <- integrate(function(t) T1_s(t)^2, 0, 1)$value |> sqrt()


T2_s <- function(x)
{
  inner1 <- integrate(function(t) (T1_s(t)/norm1_back)*(t^2), 0,1)$value
  x^2 - T1_s(x)*inner1/norm1_back - (1/3)
}

norm2_back <- integrate(function(t) T2_s(t)^2, 0, 1)$value |> sqrt()


T3_s <- function(x)
{
  inner1 <- integrate(function(t) T1_s(t)*(t^3)/norm1_back, 0,1)$value
  inner2 <- integrate(function(t) T2_s(t)*(t^3)/norm2_back, 0,1)$value
  
  x^3 - T1_s(x)*inner1/norm1_back - T2_s(x)*inner2/norm2_back - (1/4)
}

norm3_back <- integrate(function(t) T3_s(t)^2, 0, 1)$value |> sqrt()

mod_back <- function(x, beta)
{
  blen <- length(beta)
  if(blen == 1) Tvec <- T1_s(x)/norm1_back
  if(blen == 2) Tvec <- c(T1_s(x)/norm1_back, T2_s(x)/norm2_back)
  if(blen == 3) Tvec <- c(T1_s(x)/norm1_back, T2_s(x)/norm2_back,T3_s(x)/norm3_back)
  
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(back)
}

neg_loglikelihood_back <- function(theta, data)
{
  beta <- theta
  y <- sapply(data, function(t) mod_back(x = t,
                                         beta = beta))
  -sum(log(y))
}


score_back <- function(theta, data)
{
  beta <- theta
  rec <- sapply(data, function(t)
  {
    1/mod_back(x = t, beta = beta)
  })
  
  blen <- length(beta)
  if(blen == 1){
    Tmat <- sapply(data, function(t) T1_s(t)/norm1_back)
    dim(Tmat) <- c(1, length(data))
    der_beta <- Tmat%*%rec
    
  }
  if(blen == 2){
    Tmat <- rbind(sapply(data, function(t) T1_s(t)/norm1_back),
                  sapply(data, function(t) T2_s(t)/norm2_back))
    der_beta <- Tmat%*%rec
    
  }
  if(blen == 3){
    Tmat <- rbind(sapply(data, function(t) T1_s(t)/norm1_back),
                  sapply(data, function(t) T2_s(t)/norm2_back),
                  sapply(data, function(t) T3_s(t)/norm3_back))
    
    der_beta <- Tmat%*%rec
  }
  
  return(-1* as.vector(der_beta))
}