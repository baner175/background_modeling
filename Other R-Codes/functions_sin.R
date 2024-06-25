a_sig <- 6; b_sig <- 6
sig_prop <- 0.05

# make_data <- function(n = 2e3, rate_back = rate_back,
#                       a_sig = a_sig, b_sig = b_sig,
#                       bkg_prop, seed = NULL)
# {
#   library("truncdist")
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#   back <- rtrunc(n,spec="exp",rate=2,a=0,b=1)
#   sig <- rbeta(n, shape1 = a_sig, shape2 = b_sig)
#   u <- runif(n)
#   obs <- ifelse(u<bkg_prop, back,sig)
#   
#   return(list(signal = sig, background = back, observed = obs))
# }


#-------------------------------------------------------------------------------


S <- function(x, a, b)
{
  f_sig <- dbeta(x, shape1 = a, shape2 = b)
  gb <- dbeta(x, shape1 = 2, shape2 = 2)
  norm <- beta(2,2)*beta(2*a-2, 2*b-2)/((beta(a,b))^2) -1
  norm <- sqrt(norm)
  return(list('s' = f_sig/gb - 1, 's1' = (f_sig/gb-1)/norm))
}

T1 <- function(x, ...)
{
  inner <- integrate(function(t) S(t, ...)$s1*t*dbeta(t,2,2), 0, 1)$value
  x - inner*S(x, ...)$s1 - beta(3,2)/beta(2,2)
}

norm1 <- integrate(function(t) dbeta(t,2,2)*T1(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm1 <- sqrt(norm1)


T2 <- function(x, ...)
{
  inner1 <- integrate(function(t) (T1(t, ...)/norm1)*(t^2)*dbeta(t,2,2), 0,1)$value
  inner2 <- integrate(function(t) S(t, ...)$s1*(t^2)*dbeta(t,2,2), 0,1)$value
  
  x^2 - T1(x,...)*inner1/norm1 - S(x,...)$s1*inner2 - beta(4,2)/beta(2,2)
}

norm2 <- integrate(function(t) dbeta(t,2,2)*T2(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm2 <- sqrt(norm2)

T3 <- function(x, ...)
{
  inner1 <- integrate(function(t) dbeta(t,2,2)*T1(t, ...)*(t^3)/norm1, 0,1)$value
  inner2 <- integrate(function(t) dbeta(t,2,2)*T2(t, ...)*(t^3)/norm2, 0,1)$value
  inner3 <- integrate(function(t) dbeta(t,2,2)*S(t, ...)$s1*(t^3), 0,1)$value
  
  x^3 - T1(x,...)*inner1/norm1 - T2(x,...)*inner2/norm2 - S(x,...)$s1*inner3 - beta(5,2)/beta(2,2)
}

norm3 <- integrate(function(t) dbeta(t,2,2)*T3(t, a = a_sig,
                                  b = b_sig)^2, 0, 1)$value
norm3 <- sqrt(norm3)

T4 <- function(x, ...)
{
  inner1 <- integrate(function(t) dbeta(t,2,2)*T1(t, ...)*(t^4)/norm1, 0,1)$value
  inner2 <- integrate(function(t) dbeta(t,2,2)*T2(t, ...)*(t^4)/norm2, 0,1)$value
  inner3 <- integrate(function(t) dbeta(t,2,2)*T3(t, ...)*(t^4)/norm3, 0,1)$value
  inner4 <- integrate(function(t) dbeta(t,2,2)*S(t, ...)$s1*(t^4), 0,1)$value
  
  x^4 - 
    T1(x,...)*inner1/norm1 - 
    T2(x,...)*inner2/norm2 - 
    T3(x,...)*inner3/norm3 - 
    S(x,...)$s1*inner4 - beta(6,2)/beta(2,2)
}

norm4 <- integrate(function(t) dbeta(t,2,2)*T4(t, a = a_sig,
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
  gb <- dbeta(x,2,2)
  mix <- 1 + eta*S(x, ...)$s + (1-eta)*as.numeric(crossprod(beta,Tvec))
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(list(mix = gb*mix, back = gb*back))
  
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
  gb <- dbeta(data, 2,2)
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
                  sapply(data, function(t) dbeta(t,2,2)*T_basis[[j]](t, a = a_sig, b = b_sig)/norm_vec[j]))
  }
  der_beta <- (1-eta)*Tmat%*%rec
  
  beta_T_dot <- apply(Tmat, 2,
                      function(t) as.numeric(crossprod(beta, t)))
  der_eta_num <- gb*S(data, a = a_sig, b = b_sig)$s - beta_T_dot
  
  return(-1*c(
    sum(der_eta_num*rec), as.vector(der_beta)
  ))
}


mod_back_new <- function(x, beta, ...)
{
  gb <- dbeta(x,2,2)
  blen <- length(beta)
  Tvec <- c()
  for(j in 1:blen) {
    Tvec <- c(Tvec,
              T_basis[[j]](x, ...)/norm_vec[j])
  }
  back <- 1 + as.numeric(crossprod(beta,Tvec))
  return(gb*back)
}

neg_loglikelihood_back_new <- function(theta, data, a_sig, b_sig)
{
  beta <- theta
  y <- sapply(data, function(t) mod_back_new(x = t,
                                             beta = beta, a = a_sig,
                                             b = b_sig))
  -sum(log(y))
}

