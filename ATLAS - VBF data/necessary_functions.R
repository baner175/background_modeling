rm(list = ls())
library(ggplot2)
library(kdensity)
library(truncdist)
library(VGAM)
library(latex2exp)
library(Hmisc)

options(digits = 22)

mean_sig <- 125; sd_sig <- 3
eps <- 1e-3
u <- 160; l <- 110

# DEFINING SIGNAL DENSITY AND CALCULATING SIGNAL REGION

fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

# FIGURING OUT (mu_s -d, mu_s + d):
find_d <- function(d)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))

r <- sol$root

(M_lower <- mean_sig - r)
(M_upper <- mean_sig + r)

round(integrate(fs, M_lower, M_upper)$value,5) == 1-eps


# defining gb:
bkg_loc <- 91.2
bkg_scale <- 2.49/2

qb <- function(x)
{
  dtrunc(x, a = l, b= u, spec = 'cauchy', 
         location = bkg_loc, 
         scale = bkg_scale)
}

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2;

gb_test <- function(x, fs_prop = 0)
{
  fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                    spec = 'norm', a = l, b = u)
  fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                     spec = 'norm', a = l, b = u)
  qb_val <- dtrunc(x, a = l, b= u, spec = 'cauchy', 
                   location = bkg_loc, 
                   scale = bkg_scale)
  return(fs_prop*fs_val1 + fs_prop*fs_val2 + (1-2*fs_prop)*qb_val)
}

gb <- function(x) gb_test(x, fs_prop = 0.05)
gb <- Vectorize(gb)


inner_prod <- function(f1, f2, limits = c(l,u))
{
  l <- limits[1]; u <- limits[2]
  return(integrate(function(t) f1(t)*f2(t), l, u)$value)
}

# Convarting fs into S1:
norm_S <- inner_prod(function(t) (fs(t)/gb(t)-1)^2, gb, c(l,u)) |> sqrt()

S1 <- function(x)
{
  f_sig <- fs(x)
  g_b <- gb(x)
  return((f_sig/g_b -1)/norm_S)
}


construct_basis <- function(n_basis  = 2, 
                            sig_density,
                            proposal_bkg,
                            limits = c(-Inf, Inf))
{
  M <- n_basis+1
  l <- limits[1]; u <- limits[2]
  gb <- proposal_bkg; fs <- sig_density
  
  norm_S <- inner_prod(function(t) (fs(t)/gb(t)- 1)^2,
                       gb,
                       c(l,u)) |> sqrt()
  S1 <- function(x)
  {
    f_sig <- fs(x)
    g_b <- gb(x)
    return((f_sig/g_b -1)/norm_S)
  }
  S1 <- Vectorize(S1)
  
  T_basis <- list()
  T_basis[[1]] <- function(t) return(1)
  T_basis[[1]] <- Vectorize(T_basis[[1]])
  for(i in 2:M)
  {
    T_basis[[i]] <- local({
      i_val <- i
      T_inner <- c() 
      for(j in 1:(i_val-1))
      {
        T_inner[j] <- inner_prod(f1 = function(t) (t^(i_val-1))*T_basis[[j]](t),
                                 f2 = gb, limits = c(l,u))
      }
      T_inner <- c(T_inner,
                   inner_prod(f1 = function(t) (t^(i_val-1))*S1(t),
                              f2 = gb, limits = c(l,u)))
      
      Ti_unnormed <- function(t){
        fun_vals <- sapply(c(T_basis[1:(i_val-1)], S1), function(f) f(t))
        return(t^(i_val - 1) - as.numeric(crossprod(fun_vals, T_inner)))
      }
      Ti_unnormed <- Vectorize(Ti_unnormed)
      norm_Ti <- inner_prod(function(t) Ti_unnormed(t)^2,
                            gb, c(l,u)) |> sqrt()
      
      Ti <- function(t) Ti_unnormed(t)/norm_Ti
      return(Ti)
    })
  }
 return(c(T_basis[[1]], S1, T_basis[-1])) 
}
# 
# T_basis <- construct_basis(n_basis = 4,
#                            sig_density = fs, proposal_bkg = gb,
#                            limits = c(l,u))
# 
# sapply(T_basis, function(f) integrate(function(t) gb(t)*f(t)^2, l, u)$value)
# 
# curve(T_basis[[1]](x), from = l, to = u, lwd = 2, col = 'black', ylim = c(-1,3))
# curve(T_basis[[2]](x), lwd = 2, col = 'skyblue', add = TRUE)
# curve(T_basis[[3]](x), lwd = 2, col = 'red', add = TRUE)
# curve(T_basis[[4]](x), lwd = 2, col = 'blue', add = TRUE)
# curve(T_basis[[5]](x), lwd = 2, col = 'green', add = TRUE)
# curve(T_basis[[6]](x), lwd = 2, col = 'orange', add = TRUE)
