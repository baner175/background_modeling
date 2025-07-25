rm(list = ls())
library(truncdist)
library(VGAM)
library(latex2exp)
library(knitr)
library(kableExtra)
library(numDeriv)

# source('signal_search_functions.R')

### GENERIC PARAMETERS #########################################################
mean_sig <- 125; sd_sig <- 3
u <- 160; l <- 110

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
n <- length(phys_data)

file_bkg <- paste0('benchmark_toys/Cat',cat,'_',
                   scenario,'_mu0p0','.csv')
bkg_data <- read.csv(file_bkg, header = FALSE)[,1]
m <- length(bkg_data)

################################################################################



###### SPURIOUS SIGNAL ###########################################

# Here we shall estimate beta in qb for only mu = 125, and then we'll estimate
# epsilon for varying mu between 120 to 130 by fixng beta at the 
# aforementioned estimate

sp_beta_eps_mu_likelihood <- function(eps, beta, mu)
{
  qb_mass <- integrate(function(t){
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(bkg_data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  fs_i <- dtrunc(bkg_data, spec = 'norm', a = l, b = u, 
                 mean = mu, sd = sd_sig)
  fi <- (1-eps)*qb_i + eps*fs_i
  return(-sum(log(fi)))
}

sp_beta_eps_mus_likelihood <- function(par){
  eps <- par[1]; beta <- par[2]; mu <- mean_sig
  return(
    sp_beta_eps_mu_likelihood(eps, beta, mu)
  )
}

fit0 <- nlminb(start = c(0.01, 0.01),
               objective = sp_beta_eps_mus_likelihood,
               lower = c(0, -Inf),
               upper = c(1, Inf))
fun <- ifelse(fit0$convergence == 0,
              message, warning)
msg <- ifelse(fit0$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)
hess <- hessian(sp_beta_eps_mus_likelihood, x = fit0$par)

beta_hat_sp <- fit0$par[2]
beta_hat_se <- sqrt(diag(solve(hess))[2])

mu_seq <- seq(120, 130, length.out = 11)
sp_res <- sapply(mu_seq, function(t)
{
  temp_objective <- function(par){
    eps <- par[1]; beta <- par[2]; mu <- t
    return(
      sp_beta_eps_mu_likelihood(eps, beta, mu)
    )
  }
  fit <- nlminb(start = c(0.01, 0.01),
                objective = temp_objective,
                lower = c(0, -Inf),
                upper = c(1, Inf))
  hess <- hessian(temp_objective, x = fit$par)
  eps_hat <- fit$par[1]
  eps_hat_se <- sqrt(diag(solve(hess))[1])
  msg <- paste0('Evaluating at mu: ', t, '; Convergence: ', fit$convergence,'\n')
  message(msg)
  return(c(eps_hat, eps_hat_se))
})
max_idx <- which.max(sp_res[1,])
gb_sp_eps <- sp_res[1,max_idx]
gb_sp_eps_se <- sp_res[2,max_idx]

spurious_joint_likelihood_C <- function(eta, eps, beta){
  qb_mass <- integrate(function(t){
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(phys_data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  fs_i <- dtrunc(phys_data, spec = 'norm', a = l, b = u, 
                 mean = mean_sig, sd = sd_sig)
  fi <- eta*fs_i + (1-eta)*((1-eps)*qb_i + eps*fs_i)
  beta_log_lik <- dlnorm(beta, meanlog = log(beta_hat_sp), 
                         sdlog = beta_hat_se/abs(beta_hat_sp), 
                         log = TRUE)
  # eps_log_lik <- dnorm(eps, mean = gb_sp_eps,
  #                       sd = gb_sp_eps_se,
  #                       log = TRUE)
  eps_log_lik <- dlnorm(eps, meanlog = log(gb_sp_eps),
                       sdlog = gb_sp_eps_se/abs(gb_sp_eps),
                       log = TRUE)
  return(-sum(log(fi)) - beta_log_lik - eps_log_lik)
}
spurious_phys_res <- mle(minuslogl = spurious_joint_likelihood_C,
                         start = list(eta = 0.01, eps = gb_sp_eps_se, beta = beta_hat_sp),
                         method = "L-BFGS-B",
                         lower = c(0, 0, -Inf),
                         upper = c(1, 1, Inf))
fun <- ifelse(spurious_phys_res@details$convergence == 0,
              message, warning)
msg <- ifelse(spurious_phys_res@details$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

eps_phys_hat <- as.numeric(spurious_phys_res@coef['eps'])
beta_phys_hat <- as.numeric(spurious_phys_res@coef['beta'])
eta_sp <- as.numeric(spurious_phys_res@coef['eta'])
eta_sp_se <- as.numeric(summary(spurious_phys_res)@coef[,'Std. Error']['eta'])
spurious_null_res <- mle(minuslogl = spurious_joint_likelihood_C,
                         start = list(eps = gb_sp_eps_se, beta = beta_hat_sp),
                         fixed = list(eta = 0),
                         method = "L-BFGS-B",
                         lower = c(0, 0, -Inf),
                         upper = c(0, 1, Inf))
fun <- ifelse(spurious_null_res@details$convergence == 0,
              message, warning)
msg <- ifelse(spurious_null_res@details$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

ll0 <- -summary(spurious_null_res)@m2logL
ll1 <- -summary(spurious_phys_res)@m2logL
p_val_sp <- 0.5*pchisq((ll1-ll0),
                        df = 1, lower.tail = FALSE)

sp_row <- data.frame(
  'spur with Constraint', eta_sp, p_val_sp, qnorm(p_val_sp, lower.tail = FALSE) |> round(1),
  round(eta_sp*n, 2))


# qb used for safeguard
qb_sp_ <- function(x) (x^(beta_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2)
qb_sp_mass <- integrate(qb_sp_, l, u)$value
qb_sp <- function(x) (x<=u)*(x>=l)*((x^(beta_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_mass
qb_sp <- Vectorize(qb_sp)
gb_sp <- function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_mass
  
  return(eps_phys_hat*fs + (1-eps_phys_hat)*qb)
}
gb_sp <- Vectorize(gb_sp)
# curve(gb_sp, l, u)
#-------------------------------------------------------------------------------
############ OUR METHOD WITH BKG DATA ##########################################

qb_likelihood <- function(beta, data){
  qb_mass <- integrate(function(t){
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  return(-sum(log(qb_i)))
}
qb_wbkg_res <- nlminb(start = 0.01,
                      objective = qb_likelihood,
                      lower = -Inf, upper = Inf,
                      data = bkg_data)

fun <- ifelse(qb_wbkg_res$convergence == 0,
              message, warning)
msg <- ifelse(qb_wbkg_res$convergence == 0,
       'CONVERGENCE ACHIEVED!!',
       'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

beta_hat_wbkg <- qb_wbkg_res$par
qb_mass_wbkg <- integrate(function(x){
  (x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value

qb_wbkg <- function(x) ((x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg

norm_S_qb_wbkg <- integrate(function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  S_val <- (fs/qb - 1)
  return((S_val^2)*qb)
}, l, u)$value |> sqrt()

S2_phys_vec_wbkg <- sapply(phys_data, function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S_qb_wbkg^2))
})

S2_bkg_vec_wbkg <- sapply(bkg_data, function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  S_val <- (fs/qb - 1)
  return(S_val/(norm_S_qb_wbkg^2))
})

theta0_hat_wbkg <- mean(S2_phys_vec_wbkg)
delta0_hat_wbkg <- mean(S2_bkg_vec_wbkg)
eta_hat_wbkg <- (theta0_hat_wbkg - delta0_hat_wbkg)/(1-delta0_hat_wbkg)


d_log_qb_wbkg <- function(x){
  log(x) - integrate(function(t){
    qb <- ((t^(beta_hat_wbkg))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
    return(log(t)*qb)
  }, l, u)$value
}
V_hat_wbkg <- sapply(bkg_data, function(t) d_log_qb_wbkg(t)^2) |> mean()
J_hat_wbkg <- integrate(function(t) (log(t)^2)*qb_wbkg(t), l, u)$value -
  integrate(function(t) log(t)*qb_wbkg(t), l, u)$value^2

d_normS_sq_wbkg <- -integrate(function(t) d_log_qb_wbkg(t)*(fs(t)^2)/qb_wbkg(t),
                              l, u)$value
d_S2_bkg <- function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_hat_wbkg))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass_wbkg
  return(
    -((norm_S_qb_wbkg^2)*(fs/qb)*d_log_qb_wbkg(x) +
        (fs/qb - 1)*d_normS_sq_wbkg)/(norm_S_qb_wbkg^4)
  )
}

var_S2_F_hat_wbkg <- mean(S2_phys_vec_wbkg^2) - theta0_hat_wbkg^2
var_S2_Fb_hat_wbkg <- mean(S2_bkg_vec_wbkg^2) - delta0_hat_wbkg^2

test_num <- sqrt(m*n)*(eta_hat_wbkg - 0)

denom1 <- (m/((1-delta0_hat_wbkg)^2)) * var_S2_F_hat_wbkg

denom2 <- n*(((theta0_hat_wbkg-1)^2)/((1-delta0_hat_wbkg)^4)) * var_S2_Fb_hat_wbkg

denom3 <- n*(V_hat_wbkg/(J_hat_wbkg^2))*(
  (1/(1-delta0_hat_wbkg))*mean(sapply(phys_data, d_S2_bkg)) +
    ((theta0_hat_wbkg-1)/((1-delta0_hat_wbkg)^2))*mean(sapply(bkg_data, d_S2_bkg))
)^2

denom4 <- 2*(n/J_hat_wbkg)*((theta0_hat_wbkg-1)/((1-delta0_hat_wbkg)^2))*
  mean(S2_bkg_vec_wbkg*sapply(bkg_data, d_log_qb_wbkg))*
  (
    (1/(1-delta0_hat_wbkg))*mean(sapply(phys_data, d_S2_bkg)) +
      ((theta0_hat_wbkg-1)/((1-delta0_hat_wbkg)^2))*mean(sapply(bkg_data, d_S2_bkg))
  )

test_denom <- sqrt(denom1 + denom2 + denom3 + denom4)


eta_test_stat <- test_num/test_denom
p_val_eta <- pnorm(eta_test_stat, lower.tail = FALSE)

w_bkg_row <- data.frame(
  'unbiased test', eta_hat_wbkg, p_val_eta, qnorm(p_val_eta, lower.tail = FALSE) |> round(1),
  round(eta_hat_wbkg*n, 2))
#-------------------------------------------------------------------------------
############ OUR METHOD WITHOUT BKG DATA #######################################

qb_res <- nlminb(start = 0.01,
                   objective = qb_likelihood,
                   lower = -Inf, upper = Inf,
                   data = phys_data)
fun <- ifelse(qb_res$convergence == 0,
              message, warning)
msg <- ifelse(qb_res$convergence == 0,
              'CONVERGENCE ACHIEVED!!',
              'XXX !!! FAILED TO CONVERGE !!! XXX')
fun(msg)

beta_hat <- qb_res$par

qb_mass <- integrate(function(x){
  (x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2)
}, l, u)$value
qb <- function(x) ((x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass

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
    qb <- ((x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    return(((fs/gb - 1)^2)*gb)
  }, l, u)$value |> sqrt()
  
  S2_phys_vec <- sapply(phys_data, function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
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
  
  theta0_hat <- mean(S2_phys_vec)
  
  lgx_qb_int <- integrate(function(t){
    qb <- ((t^(beta_hat))/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
    return(log(t)*qb)
  }, l, u)$value
  
  d_norm_S_sq <- -(1-2*lambda)*integrate(function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- log(x) - lgx_qb_int
    return(((fs/gb)^2)*qb*d_log_qb)
  }, l, u)$value
  
  d_S2 <- function(x){
    fs <- dtrunc(x, spec = 'norm', a = l, b = u,
                 mean = mean_sig, sd = sd_sig)
    qb <- ((x^(beta_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mass
    fs_val1 <- dtrunc(x, mean = mean1_in_gb, sd = sd_in_gb,
                      a = l, b = u,
                      spec = 'norm')
    fs_val2 <-  dtrunc(x, mean = mean2_in_gb, sd = sd_in_gb,
                       a = l, b = u,
                       spec = 'norm')
    gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*qb
    d_log_qb <- log(x) - lgx_qb_int
    
    num1 <- (norm_S^2)*((fs*qb)/(gb^2))*(1-2*lambda)*d_log_qb
    num2 <- (fs/gb - 1)*d_norm_S_sq
    return(
      -(num1+num2)/(norm_S^4)
    )
  }
  var_S2_F_hat <- mean(S2_phys_vec^2) - theta0_hat^2
  
  V_hat <- sapply(phys_data, function(x){
    (log(x) - lgx_qb_int)^2
  }) |> mean()
  J_hat <- integrate(function(t) (log(t)^2)*qb(t), l, u)$value - 
    integrate(function(t) log(t)*qb(t), l, u)$value^2
  
  sig_theta0_hat <- sqrt(
    var_S2_F_hat + (V_hat/(J_hat^2))*mean(sapply(phys_data, d_S2))^2 + 
      (2/J_hat)*mean(sapply(phys_data, d_S2))*
      mean(S2_phys_vec*sapply(phys_data, function(x){log(x) - lgx_qb_int}))
  )
  
  theta0_stat <- sqrt(n)*(theta0_hat-0)/sig_theta0_hat
  p_val <- pnorm(theta0_stat, lower.tail = FALSE)
  S_hat <- n*theta0_hat
  
  return(c(theta0_hat,
           p_val,
           qnorm(p_val, lower.tail = FALSE) |> round(1),
           round(S_hat, 2)))
}

lambda_max <- 0.05 # change lambda_max here
lambda_seq <- seq(0, lambda_max, length.out = 5)

res_sig_search <- sapply(lambda_seq, signal_search)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
res_sig_search[,1] <- as.character(res_sig_search[,1])

################################################################################

######## Creating results table ################################################

colnames(res_sig_search) <- colnames(sp_row) <- colnames(w_bkg_row) <- 
  c("$\\lambda$ / Test", "$\\hat{\\eta}$", "$p$-value", '$\\sigma$-signif',
    '$\\hat{S}$')

res_sig_search <- rbind(res_sig_search, w_bkg_row, sp_row)

caption <- paste0("Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2, '<br> Bkg: BW X Power0')
kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center')

################################################################################

######### Plotting Densities ###################################################

hist(phys_data, probability = TRUE, breaks = 50,
     col = 'white', 
     main = TeX('Densities $g^{spur}_b(x;\\hat{\\epsilon},\\hat{\\beta})$ and $q_b(x; \\hat{beta}_{bkg})$'))
curve(gb_sp, l, u, add = TRUE,
      col = 'blue', lwd = 2)
curve(qb_sp, l, u, add = TRUE,
      col = 'purple', lwd = 2, lty = 2)
legend('topright', col = c('blue', 'purple'),
       lty = 1:2, bty = 'n', lwd = 2.2,
       legend=c(TeX('$g^{spur}_b(x; \\hat{\\epsilon},\\hat{\\beta})$'),
                TeX('$q_b(x; \\hat{beta}_{bkg})$')),
       cex = 1,
       y.intersp = 1)

picture_l <- M_lower - 3; picture_u <- M_upper + 3

mycols <- c('brown', 'skyblue' ,'red', 'orange', 'purple')
palette(mycols)
my_lty = 1:5

hist(phys_data, probability = TRUE, breaks = 50,
     xlim = c(picture_l, picture_u), col = 'white',
     main = 'Sensitivity analysis without bkg data')

for(j in 1:length(lambda_seq))
{
  curve((dtrunc(x, spec = 'norm', a = l, b = u,
                mean = mean1_in_gb, sd = sd_in_gb) + 
           dtrunc(x, spec = 'norm', a = l, b = u,
                  mean = mean2_in_gb, sd = sd_in_gb))*lambda_seq[j] +
          (1-2*lambda_seq[j])*qb(x),
        l, u, add = TRUE, lwd = 2.2,
        col = ggplot2::alpha(palette()[j], alpha = 0.6),
        lty = my_lty[j])
}

legend('topright', col = mycols,
       lty = 1:5, bty = 'n', lwd =2.2,
       legend=TeX(sprintf(r'($g_b(\lambda = %f)$)', lambda_seq)),
       cex = 1,
       y.intersp = 1)
