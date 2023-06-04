a_back <- 1; b_back <- 5
a_sig <- 16; b_sig <- 2
sig_prop <- 0.05

make_data <- function(n = 2e3, a_back = a_back, 
                      b_back = b_back,
                      a_sig = a_sig, b_sig = b_sig,
                      bkg_prop, seed = NULL)
{
  if(!is.null(seed)){
    set.seed(seed)
  }
  back <- rbeta(n, shape1 = a_back, shape2 = b_back)
  sig <- rbeta(n, shape1 = a_sig, shape2 = b_sig)
  u <- runif(n)
  obs <- ifelse(u<bkg_prop, back,sig)
  
  return(list(signal = sig, background = back, observed = obs))
}


#-------------------------------------------------------------------------------


S <- function(x, a, b)
{
  f_sig <- dbeta(x, shape1 = a, shape2 = b)
  norm <- beta(2*a-1, 2*b-1)/((beta(a,b))^2) -1
  norm <- sqrt(norm)
  return(list('s' = f_sig - 1, 's1' = (f_sig-1)/norm))
}

T1 <- function(x, ...)
{
  inner <- integrate(function(t) S(t, ...)$s1*t, 0, 1)$value
  x - inner*S(x, ...)$s1 - 0.5
}

norm1 <- integrate(function(t) T1(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm1 <- sqrt(norm1)


T2 <- function(x, ...)
{
  inner1 <- integrate(function(t) (T1(t, ...)/norm1)*(t^2), 0,1)$value
  inner2 <- integrate(function(t) S(t, ...)$s1*(t^2), 0,1)$value
  
  x^2 - T1(x,...)*inner1/norm1 - S(x,...)$s1*inner2 - (1/3)
}

norm2 <- integrate(function(t) T2(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm2 <- sqrt(norm2)

T3 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^3)/norm1, 0,1)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^3)/norm2, 0,1)$value
  inner3 <- integrate(function(t) S(t, ...)$s1*(t^3), 0,1)$value
  
  x^3 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - S(x,...)$s1*inner3 - (1/4)
}

norm3 <- integrate(function(t) T3(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm3 <- sqrt(norm3)

T4 <- function(x, ...)
{
  inner1 <- integrate(function(t) T1(t, ...)*(t^4)/norm1, 0,1)$value
  inner2 <- integrate(function(t) T2(t, ...)*(t^4)/norm2, 0,1)$value
  inner3 <- integrate(function(t) T3(t, ...)*(t^4)/norm3, 0,1)$value
  inner4 <- integrate(function(t) S(t, ...)$s1*(t^4), 0,1)$value
  
  x^4 - 
    T1(x,...)*inner1/norm1 - 
    T2(x,...)*inner2/norm2 - 
    T3(x,...)*inner3/norm3 - 
    S(x,...)$s1*inner4 - (1/5)
}

norm4 <- integrate(function(t) T4(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
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

neg_loglikelihood <- function(theta, data, a_sig, b_sig)
{
  eta = theta[1]
  beta = theta[-1]
  y <- sapply(data, function(t) mod(x = t, eta = eta,
                                    beta = beta, a = a_sig, 
                                    b = b_sig)$mix)
  -sum(log(y))
}

score <- function(theta, data, a_sig, b_sig)
{
  eta = theta[1]
  beta = theta[-1]
  
  rec <- sapply(data, function(t)
  {
    1/mod(x = t,  eta = eta,
          beta = beta, a = a_sig,
          b = b_sig)$mix
  })
  
  blen <- length(beta)
  Tmat <- c()
  for(j in 1:blen)
  {
    Tmat <- rbind(Tmat,
                  sapply(data, function(t) T_basis[[j]](t, a = a_sig, b = b_sig)/norm_vec[j]))
  }
  der_beta <- (1-eta)*Tmat%*%rec
  
  der_eta_num <- sapply(data, function(t)
  {
    index <- which(data == t)
    Tvec <- Tmat[,index]
    return(S(t, a = a_sig, b = b_sig)$s -
             as.numeric(crossprod(beta, Tvec)))
  })
  
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

neg_loglikelihood_back_new <- function(theta, data, a_sig, b_sig)
{
  beta <- theta
  y <- sapply(data, function(t) mod_back_new(x = t,
                                             beta = beta, a = a_sig,
                                             b = b_sig))
  -sum(log(y))
}

