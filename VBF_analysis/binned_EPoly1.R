rm(list = ls())
library(truncdist)
library(VGAM)
library(latex2exp)
library(knitr)
library(kableExtra)
library(numDeriv)
library(nloptr)

# source('signal_search_functions.R')

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110
k <- 5e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

# Defining signal density and calculating signal region
fs <- function(x) dtrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)
Fs <- function(x) ptrunc(x, spec = 'norm', a = l, b = u, mean = mean_sig, sd = sd_sig)

# Defining BW parameters
bkg_loc <- 91.2
bkg_scale <- 2.49/2

### LOADING DATA ###############################################################
cat <- 0
scenario <- 'HLHC'
mu <- 1

mu_part_1 <- floor(mu); mu_part_2 <- (10*mu)%%10

file_name <- paste0('benchmark_toys/Cat',cat,'_',
                    scenario,'_mu',mu_part_1,'p',mu_part_2,'.csv')

phys_data <- read.csv(file_name, header = FALSE)[,1]
N <- length(phys_data)
ni <- sapply(1:k, function(i){
  sum((phys_data>bin_ends[i])&(phys_data<=bin_ends[i+1]))
})

file_bkg <- paste0('benchmark_toys/Cat',cat,'_',
                   scenario,'_mu0p0','.csv')
bkg_data <- read.csv(file_bkg, header = FALSE)[,1]
M <- length(bkg_data)
mi <- sapply(1:k, function(i){
  sum((bkg_data>bin_ends[i])&(bkg_data<=bin_ends[i+1]))
})

################################################################################
### ESTIMATING BETA WITH BACKGROUND ONLY DATA ##################################

qb_likelihood <- function(beta, bin_counts){
  qb_mass <- integrate(function(t){
    (exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(1:k, function(i){
    integrate(function(t){
      ((exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(bin_counts*log(qb_i)))
}
qb_wbkg_res <- nlminb(start = 0.01,
                      objective = qb_likelihood,
                      lower = -Inf, upper = Inf,
                      bin_counts = mi)

fun <- ifelse(qb_wbkg_res$convergence == 0,
              message, warning)
msg <- ifelse(qb_wbkg_res$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

beta_hat_wbkg <- qb_wbkg_res$par
beta_hat_wbkg_se <- (1/sqrt(hessian(qb_likelihood, 
                                    x = beta_hat_wbkg, 
                                    bin_counts = mi)))
################################################################################
#-------------------------------------------------------------------------------
###### SPURIOUS SIGNAL MAX SPUR ##########################################

# Here we shall proceed with the estimated beta using only a qb-likelihood
# and then we'll estimate epsilon for varying mu between 120 to 130 
# by fixng beta at the aforementioned estimate

sp_eps_mu_likelihood <- function(eps, mu)
{
  beta <- beta_hat_wbkg
  qb_mass <- integrate(function(t){
    (exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(1:k, function(i){
    integrate(function(t){
      ((exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fs_i <- sapply(1:k, function(i){
    integrate(function(t){
      dtrunc(t, spec = 'norm', a = l, b = u,
             mean = mu, sd = sd_sig)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fi <- (1-eps)*qb_i + eps*fs_i
  return(-sum(mi*log(fi)))
}

mu_seq <- seq(120, 130, length.out = 11)
sp_res <- sapply(mu_seq, function(t)
{
  temp_objective <- function(par){
    eps <- par; mu <- t
    return(
      sp_eps_mu_likelihood(eps, mu)
    )
  }
  fit <- nlminb(start = 0.01,
                objective = temp_objective,
                lower = 0,
                upper = 1)
  hess <- hessian(temp_objective, x = fit$par)
  eps_hat <- fit$par[1]
  eps_hat_se <- sqrt(diag(solve(hess))[1])
  msg <- paste0('Evaluating at mu: ', t, '; Convergence: ', fit$convergence,'\n')
  message(msg)
  return(c(eps_hat, eps_hat_se))
})
max_idx <- which.max(sp_res[1,])
eps_sp_max_hat <- sp_res[1,max_idx]
eps_sp_max_hat_se <- sp_res[2,max_idx]

spurious_max_joint_likelihood_C <- function(par){
  eta <- par[1]; eps <- par[2]; beta <- par[3]
  qb_mass <- integrate(function(t){
    (exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  
  qb_i <- sapply(1:k, function(i){
    integrate(function(t){
      ((exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fs_i <- sapply(1:k, function(i){
    integrate(function(t){
      dtrunc(t, spec = 'norm', a = l, b = u,
             mean = mean_sig, sd = sd_sig)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fi <- eta*fs_i + (1-eta)*((1-eps)*qb_i + eps*fs_i)
  beta_log_lik <- dlnorm(beta, meanlog = log(beta_hat_wbkg), 
                         sdlog = beta_hat_wbkg_se/abs(beta_hat_wbkg), 
                         log = TRUE)
  # eps_log_lik <- dnorm(eps, mean = eps_sp_max_hat,
  #                       sd = eps_sp_max_hat_se,
  #                       log = TRUE)
  eps_log_lik <- dlnorm(eps, meanlog = log(eps_sp_max_hat),
                        sdlog = eps_sp_max_hat_se/abs(eps_sp_max_hat),
                        log = TRUE)
  return(-sum(ni*log(fi)) - beta_log_lik - eps_log_lik)
}

opts <- list("algorithm" = "NLOPT_LN_COBYLA",
             "xtol_rel" = 1.0e-8,
             "maxeval" = 1e4)

spurious_max_phys_res <- nloptr(x0 = c(0, eps_sp_max_hat, beta_hat_wbkg),
                                eval_f = spurious_max_joint_likelihood_C,
                                opts = opts,
                                lb = c(0, 0, -Inf),
                                ub = c(1, 1, Inf))

spurious_max_joint_likelihood_null <- function(par)
{
  par_new <- c(0, par)
  return(spurious_max_joint_likelihood_C(par_new))
}

spurious_max_null_res <- nloptr(x0 = c(eps_sp_max_hat, beta_hat_wbkg),
                                eval_f = spurious_max_joint_likelihood_null,
                                opts = opts,
                                lb = c(0, -Inf),
                                ub = c(1, Inf))

eps_sp_max_phys_hat <- spurious_max_phys_res$solution[2]
beta_sp_max_phys_hat <- spurious_max_phys_res$solution[3]
eta_sp_max_hat <- spurious_max_phys_res$solution[1]
ll0 <- -spurious_max_null_res$objective
ll1 <- -spurious_max_phys_res$objective
p_val_sp_max <- 0.5*pchisq(2*(ll1-ll0),
                           df = 1, lower.tail = FALSE)

sp_max_row <- data.frame(
  'Max spur with Constraint', 
  eta_sp_max_hat, 
  p_val_sp_max,
  max(qnorm(p_val_sp_max, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_sp_max_hat*N, 2), 0))

qb_sp_max_ <- function(x) (exp(x*beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2)
qb_sp_max_mass <- integrate(qb_sp_max_, l, u)$value
qb_sp_max <- function(x) (x<=u)*(x>=l)*((exp(x*beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_max_mass
qb_sp_max <- Vectorize(qb_sp_max)
gb_sp_max <- function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((exp(x*beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_max_mass
  
  return(eps_sp_max_phys_hat*fs + (1-eps_sp_max_phys_hat)*qb)
}
gb_sp_max <- Vectorize(gb_sp_max)

#-------------------------------------------------------------------------------
############ OUR METHOD WITH BKG DATA ##########################################

qb_mass_wbkg <- integrate(function(x){
  (exp(x*beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_wbkg <- function(x) ((exp(x*beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg

norm_S_qb_wbkg <- integrate(function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((exp(x*beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  S_val <- (fs/qb - 1)
  return((S_val^2)*qb)
}, l, u)$value |> sqrt()

S2_vec_wbkg <- sapply(xi, function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((exp(x*beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S_qb_wbkg^2))
})

theta0_hat_wbkg <- sum(ni*S2_vec_wbkg)/N
delta0_hat_wbkg <- sum(mi*S2_vec_wbkg)/M
eta_hat_wbkg <- (theta0_hat_wbkg - delta0_hat_wbkg)/(1-delta0_hat_wbkg)
c_hat <- N/k; cb_hat <- M/k

d_log_qb_wbkg <- function(x){
  x - integrate(function(t){
    qb <- ((exp(t*beta_hat_wbkg))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
    return(t*qb)
  }, l, u)$value
}
d_log_qb_wbkg_xi <- sapply(xi, d_log_qb_wbkg)
d2_int1_wbkg <- integrate(function(t){
  qb <- ((exp(t*beta_hat_wbkg))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  return((t^2)*qb)
}, l, u)$value
d2_int2_wbkg <- integrate(function(t){
  qb <- ((exp(t*beta_hat_wbkg))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  return(t*qb)
}, l, u)$value
d2_log_qb_wbkg <- -(d2_int1_wbkg - d2_int2_wbkg^2)
V_hat_wbkg <- sum((d_log_qb_wbkg_xi^2)*mi)/M
J_hat_wbkg <- -sum(d2_log_qb_wbkg*mi)/k
d_normS_sq_wbkg <- -integrate(function(t) d_log_qb_wbkg(t)*(fs(t)^2)/qb_wbkg(t),
                              l, u)$value
d_S2_vec_wbkg <- sapply(xi, function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((exp(x*beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  return(
    -((norm_S_qb_wbkg^2)*(fs/qb)*d_log_qb_wbkg(x) +
        (fs/qb - 1)*d_normS_sq_wbkg)/(norm_S_qb_wbkg^4)
  )
})
d_theta0_hat_wbkg <- sum(d_S2_vec_wbkg*ni)/N
d_delta0_hat_wbkg <- sum(d_S2_vec_wbkg*mi)/M

var_S2_F_hat_wbkg <- sum((S2_vec_wbkg^2)*ni)/N - theta0_hat_wbkg^2
var_S2_Fb_hat_wbkg <- sum((S2_vec_wbkg^2)*mi)/M - delta0_hat_wbkg^2

d_theta0_T_wbkg <- 1/(1-delta0_hat_wbkg)
d_delta0_T_wbkg <- (theta0_hat_wbkg-1)/((1-delta0_hat_wbkg)^2)
cov_term_wbkg <- (sum(mi*S2_vec_wbkg*d_log_qb_wbkg_xi)/M)


test_num <- sqrt(M*N)*(eta_hat_wbkg - 0)

denom1 <- (M/((1-delta0_hat_wbkg)^2)) * var_S2_F_hat_wbkg

denom2 <- N*(((theta0_hat_wbkg-1)^2)/((1-delta0_hat_wbkg)^4)) * var_S2_Fb_hat_wbkg

denom3 <- N*(V_hat_wbkg/(J_hat_wbkg^2))*(cb_hat^2)*((d_theta0_T_wbkg*d_theta0_hat_wbkg + d_delta0_T_wbkg*d_delta0_hat_wbkg)^2)

denom4 <- (2*N*cb_hat/J_hat_wbkg)*d_delta0_T_wbkg*(d_theta0_T_wbkg*d_theta0_hat_wbkg + 
                                                     d_delta0_T_wbkg*d_delta0_hat_wbkg)*cov_term_wbkg


test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)

eta_test_stat <- test_num/test_denom
p_val_eta <- pnorm(eta_test_stat, lower.tail = FALSE)

unbiased_row <- data.frame(
  'unbiased test',
  eta_hat_wbkg, 
  p_val_eta, 
  max(qnorm(p_val_eta, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_hat_wbkg*N, 2), 0))

#-------------------------------------------------------------------------------
###### REGULAR MLE ###########################################

mle_likelihood_C <- function(par){
  eta <- par[1]; beta <- par[2]
  qb_mass <- integrate(function(t){
    (exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(1:k, function(i){
    integrate(function(t){
      ((exp(t*beta))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fs_i <- sapply(1:k, function(i){
    integrate(function(t){
      dtrunc(t, spec = 'norm', a = l, b = u, 
             mean = mean_sig, sd = sd_sig)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  fi <- eta*fs_i + (1-eta)*qb_i
  beta_log_lik <- dlnorm(beta, meanlog = log(beta_hat_wbkg), 
                         sdlog = beta_hat_wbkg_se/abs(beta_hat_wbkg), 
                         log = TRUE)
  return(-sum(ni*log(fi)) - beta_log_lik)
}

mle_phys_res <- nloptr(x0 = c(0.01, beta_hat_wbkg),
                       eval_f = mle_likelihood_C,
                       opts = opts,
                       lb = c(0, -Inf),
                       ub = c(1, Inf))

eta_mle_hat <- mle_phys_res$solution[1]
beta_mle_hat <- mle_phys_res$solution[2]


mle_likelihood_null <- function(par)
{
  par_new <- c(0, par)
  return(mle_likelihood_C(par_new))
}

mle_null_res <- nloptr(x0 = beta_hat_wbkg,
                       eval_f = mle_likelihood_null,
                       opts = opts,
                       lb = -Inf,
                       ub = Inf)

ll0 <- -mle_null_res$objective
ll1 <- -mle_phys_res$objective
p_val_mle <- 0.5*pchisq(2*(ll1-ll0),
                        df = 1, lower.tail = FALSE)

mle_row <- data.frame(
  'Regular MLE with Constraint', 
  eta_mle_hat, 
  p_val_mle, 
  max(qnorm(p_val_mle, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_mle_hat*N, 2), 0))

# qb used for safeguard
qb_mle_ <- function(x) (exp(x*beta_mle_hat))/((x-bkg_loc)^2 + bkg_scale^2)
qb_mle_mass <- integrate(qb_mle_, l, u)$value
qb_mle <- function(x) (x<=u)*(x>=l)*((exp(x*beta_mle_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mle_mass
qb_mle <- Vectorize(qb_mle)
#-------------------------------------------------------------------------------
############ OUR METHOD WITHOUT BKG DATA #######################################

qb_res <- nlminb(start = 0.01,
                 objective = qb_likelihood,
                 lower = -Inf, upper = Inf,
                 bin_counts = ni)
fun <- ifelse(qb_res$convergence == 0,
              message, warning)
msg <- ifelse(qb_res$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

beta_hat <- qb_res$par

qb_mass <- integrate(function(x){
  (exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value
qb <- function(x) ((exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass

# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d, eps = 1e-3)
{
  pl <- Fs(mean_sig-d)
  pu <- Fs(mean_sig+d)
  return(pu-pl-1+eps)
}
sol <- uniroot(find_d, lower = 0, upper = min(mean_sig - l,u - mean_sig))
r <- sol$root
M_lower <- mean_sig - r; M_upper <- mean_sig + r

# Defining gb:
mean1_in_gb <- (M_lower + mean_sig)/2; sd_in_gb <- 2*sd_sig
mean2_in_gb <- (M_upper + mean_sig)/2

signal_search <- function(lambda){
  norm_S <- integrate(function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    return(((fs/gb - 1)^2)*gb)
  }, l, u)$value |> sqrt()
  
  S2_vec <- sapply(xi, function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    S_val <- (fs/gb - 1)
    return(S_val/(norm_S^2))
  })
  d_log_qb <- function(x){
    x - integrate(function(t){
      qb <- ((exp(t*beta_hat))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
      t*qb
    }, l, u)$value
  }
  d_log_qb_xi <- sapply(xi, d_log_qb)
  d_norm_S_sq <- -(1-2*lambda)*integrate(function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- d_log_qb(x)
    return(((fs/gb)^2)*qb*d_log_qb)
  }, l, u)$value
  
  theta0_hat <- sum(ni*S2_vec)/N
  
  d2_int1 <- integrate(function(t){
    qb <- ((exp(t*beta_hat))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    (t^2)*qb
  }, l, u)$value
  d2_int2 <- integrate(function(t){
    qb <- ((exp(t*beta_hat))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    t*qb
  }, l, u)$value
  d2_log_qb <- -(d2_int1 - (d2_int2^2))
  d_S2_vec <- sapply(xi, function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((exp(x*beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- d_log_qb(x)
    num1 <- (norm_S^2)*((fs*qb)/(gb^2))*(1-2*lambda)*d_log_qb
    num2 <- (fs/gb - 1)*d_norm_S_sq
    return(
      -(num1+num2)/(norm_S^4)
    )
  })
  
  d_theta0_hat <- sum(d_S2_vec*ni)/N
  
  var_S2_F_hat <- sum((S2_vec^2)*ni)/N - theta0_hat^2
  
  V_hat <- sum(ni*(d_log_qb_xi^2))/N
  J_hat <- (-1/k)*sum(ni*d2_log_qb)
  
  sig_theta0_hat <- sqrt(
    var_S2_F_hat + 
      (c_hat^2)*(V_hat/(J_hat^2))*(d_theta0_hat^2) + 
      (2*c_hat/J_hat)*d_theta0_hat*
      (sum(ni*S2_vec*d_log_qb_xi)/N)
  )
  
  theta0_stat <- sqrt(N)*(theta0_hat-0)/sig_theta0_hat
  p_val <- pnorm(theta0_stat, lower.tail = FALSE)
  S_hat <- N*theta0_hat
  
  return(c(theta0_hat,
           p_val,
           max(qnorm(p_val, lower.tail = FALSE) |> round(1), 0),
           max(round(S_hat, 2),0)))
}

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 5)

res_sig_search <- sapply(lambda_seq, signal_search)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
res_sig_search[,1] <- as.character(res_sig_search[,1])

################################################################################

######## Creating results table ################################################

colnames(res_sig_search) <- colnames(sp_max_row) <- colnames(mle_row) <- colnames(unbiased_row) <- 
  c("$\\lambda$ / Test", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif',
    '$\\hat{S}$')

res_sig_search <- rbind(res_sig_search,
                        unbiased_row,
                        mle_row, 
                        sp_max_row)

caption <- paste0("Binned (", k," bins) Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2, '<br> Bkg: BW X EPoly1')

kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center') %>%
  row_spec(c(6,7,8), extra_css = "border-top: 2px solid black;")

################################################################################
######### Plotting Densities ###################################################

k_graph <- 100
bin_ends_graph <- seq(l, u, length.out = k_graph + 1)
xi_graph <- (bin_ends_graph[1:k_graph] + bin_ends_graph[2:(k_graph+1)])/2
ni_graph <- sapply(1:k_graph, function(i){
  sum((phys_data>bin_ends_graph[i])&(phys_data<=bin_ends_graph[i+1]))
})
mi_graph <- sapply(1:k_graph, function(i){
  sum((bkg_data>bin_ends_graph[i])&(bkg_data<=bin_ends_graph[i+1]))
})


plot(x = xi_graph, y = mi_graph, pch = 16, col = 'grey',
     main = TeX('Densities $g^{spur-max}_b(x;\\hat{\\epsilon},\\hat{\\beta})$ and $q_b^{MLE}(x; \\hat{beta}_{bkg})$ - on bkg'))

curve(M*gb_sp_max(x)*(u-l)/(k_graph+1), l, u, add = TRUE,
      col = ggplot2::alpha('purple', 0.6), lwd = 2, lty = 2)
curve(M*qb_mle(x)*(u-l)/(k_graph+1), l, u, add = TRUE,
      col = ggplot2::alpha('black', 0.6), lwd = 2, lty = 1)
legend('topright',
       col = c('purple', 'black'),
       lty = 2:1,
       bty = 'n',
       lwd = 2.2,
       legend = c(
         TeX(sprintf(r'($g^{spur-max}_b(x; \hat{\epsilon} = %.4f, \hat{\beta})$)', eps_sp_max_hat)),
         TeX(r'($q_b^{MLE}(x; \hat{\beta}_{bkg})$)')
       ),
       cex = 1,
       y.intersp = 1)

# picture_l <- M_lower - 3; picture_u <- M_upper + 3
picture_l <- 115; picture_u <- 135
mycols <- c('black', 'skyblue' ,'red', 'orange', 'brown')
palette(mycols)
my_lty = 1:5
plot(x = xi_graph, y = ni_graph, pch = 16, col = 'grey', main = '')


for(j in 1:length(lambda_seq))
{
  curve(((dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean1_in_gb, sd = sd_in_gb) +
            dtrunc(x, spec = 'norm', a = l, b = u,
                   mean = mean2_in_gb, sd = sd_in_gb))*lambda_seq[j] +
           (1-2*lambda_seq[j])*qb(x))*N*(u-l)/(k_graph+1),
        l, u, add = TRUE, lwd = 2.2,
        col = ggplot2::alpha(palette()[j], alpha = 0.6),
        lty = my_lty[j])
}

curve(N*gb_sp_max(x)*(u-l)/(k_graph+1), l, u, add = TRUE,
      col = 'purple', lwd = 2, lty = 6)
title(paste0('Sensitivity analysis + Spurious signal background:\n\n',
             'DATA: Cat',cat,'_',scenario,
             '_mu',mu_part_1,'p',mu_part_2, '; BKG: BW X EPoly1'))
legend('topright',
       col = c(mycols, 'purple'),
       lty = c(1:5, 6),
       bty = 'n',
       lwd = 2.2,
       legend = c(
         TeX(sprintf(r'($g_b(\lambda = %f)$)', lambda_seq)),
         TeX(sprintf(r'($g^{spur-max}_b(x; \hat{\epsilon} = %.4f, \hat{\beta})$)', eps_sp_max_hat))
       ),
       cex = 1,
       y.intersp = 1)
