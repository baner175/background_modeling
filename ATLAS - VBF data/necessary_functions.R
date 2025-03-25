generate_unbinned_data <- function(bin_counts, bin_ends, seed = NULL)
{
  set.seed(seed)
  bins <- bin_ends; ni <- bin_counts
  obs <- c()
  for(i in 1:(length(bins))-1)
  {
    obs <- c(obs, runif(ni[i], bins[i], bins[i+1]))
  }
  return(obs)
}

find_signal_region <- function(signal_CDF, signal_location, signal_mass)
{
  Fs <- signal_CDF
  find_d <- function(d)
  {
    pl <- Fs(signal_location-d)
    pu <- Fs(signal_location+d)
    return(pu-pl-signal_mass)
  }
  
  sol <- uniroot(find_d, lower = 0, upper = min(signal_location - l,
                                                u - signal_location))
  
  r <- sol$root
  
  return(signal_location +c(-1,1)* r)
}

inner_prod <- function(f1, f2, limits = c(l,u))
{
  l <- limits[1]; u <- limits[2]
  return(integrate(function(t) f1(t)*f2(t), l, u)$value)
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
