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
### ESTIMATING BETA WITH BACKGROUND ONLY DATA ##################################

qb_likelihood <- function(beta, data){
  qb_mass <- integrate(function(t){
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  return(-mean(log(qb_i)))
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
beta_hat_wbkg_se <- (1/sqrt(hessian(qb_likelihood, 
                                    x = beta_hat_wbkg, 
                                    data = bkg_data)))/sqrt(m)
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
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(bkg_data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  fs_i <- dtrunc(bkg_data, spec = 'norm', a = l, b = u, 
                 mean = mu, sd = sd_sig)
  fi <- (1-eps)*qb_i + eps*fs_i
  return(-mean(log(fi)))
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
  eps_hat_se <- sqrt(diag(solve(hess))[1])/sqrt(m)
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
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(phys_data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  fs_i <- dtrunc(phys_data, spec = 'norm', a = l, b = u, 
                 mean = mean_sig, sd = sd_sig)
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
  return(-sum(log(fi)) - beta_log_lik - eps_log_lik)
}

opts <- list("algorithm" = "NLOPT_LN_COBYLA",
             "xtol_rel" = 1.0e-8,
             "maxeval" = 1e4)

spurious_max_phys_res <- nloptr(x0 = c(0.01, eps_sp_max_hat, beta_hat_wbkg),
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
  max(round(eta_sp_max_hat*n, 2), 0))

qb_sp_max_ <- function(x) (x^(beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2)
qb_sp_max_mass <- integrate(qb_sp_max_, l, u)$value
qb_sp_max <- function(x) (x<=u)*(x>=l)*((x^(beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_max_mass
qb_sp_max <- Vectorize(qb_sp_max)
gb_sp_max <- function(x){
  fs <- dtrunc(x, spec = 'norm', a = l, b = u,
               mean = mean_sig, sd = sd_sig)
  qb <- (x<=u)*(x>=l)*((x^(beta_sp_max_phys_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_sp_max_mass
  
  return(eps_sp_max_phys_hat*fs + (1-eps_sp_max_phys_hat)*qb)
}
gb_sp_max <- Vectorize(gb_sp_max)

#-------------------------------------------------------------------------------
############ OUR METHOD WITH BKG DATA ##########################################

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

unbiased_row <- data.frame(
  'unbiased test',
  eta_hat_wbkg, 
  p_val_eta, 
  max(qnorm(p_val_eta, lower.tail = FALSE) |> round(1), 0),
  max(round(eta_hat_wbkg*n, 2), 0))

#-------------------------------------------------------------------------------
###### REGULAR MLE ###########################################

mle_likelihood_C <- function(par){
  eta <- par[1]; beta <- par[2]
  qb_mass <- integrate(function(t){
    (t^beta)/((t-bkg_loc)^2 + bkg_scale^2)
  }, l, u)$value
  qb_i <- sapply(phys_data, function(t){
    ((t^beta)/((t-bkg_loc)^2 + bkg_scale^2))/qb_mass
  })
  fs_i <- dtrunc(phys_data, spec = 'norm', a = l, b = u, 
                 mean = mean_sig, sd = sd_sig)
  fi <- eta*fs_i + (1-eta)*qb_i
  beta_log_lik <- dlnorm(beta, meanlog = log(beta_hat_wbkg), 
                         sdlog = beta_hat_wbkg_se/abs(beta_hat_wbkg), 
                         log = TRUE)
  return(-sum(log(fi)) - beta_log_lik)
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
  max(round(eta_mle_hat*n, 2), 0))

# qb used for safeguard
qb_mle_ <- function(x) (x^(beta_mle_hat))/((x-bkg_loc)^2 + bkg_scale^2)
qb_mle_mass <- integrate(qb_mle_, l, u)$value
qb_mle <- function(x) (x<=u)*(x>=l)*((x^(beta_mle_hat))/((x-bkg_loc)^2 + bkg_scale^2))/qb_mle_mass
qb_mle <- Vectorize(qb_mle)
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

caption <- paste0("Unbinned Signal Search Results: ", 
                  'Cat',cat,'_',scenario,
                  '_mu',mu_part_1,'p',mu_part_2, '<br> Bkg: BW X Power0')

kable(res_sig_search, "html", digits = 4,
      booktabs = TRUE, escape = FALSE,
      caption = caption) %>%
  kable_styling(full_width = FALSE, position = 'center') %>%
  row_spec(c(6,7,8), extra_css = "border-top: 2px solid black;")

################################################################################
######### Plotting Densities ###################################################

hist(bkg_data, probability = TRUE, breaks = 50,
     col = 'white', 
     main = TeX('Densities $g^{spur-max}_b(x;\\hat{\\epsilon},\\hat{\\beta})$ and $q_b^{MLE}(x; \\hat{beta}_{bkg})$ - on bkg'))
curve(gb_sp_max, l, u, add = TRUE,
      col = ggplot2::alpha('purple', 0.6), lwd = 2, lty = 2)
curve(qb_mle, l, u, add = TRUE,
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
hist(phys_data, probability = TRUE, breaks = 50,
     xlim = c(picture_l, picture_u), col = 'white',
     main = '')

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

curve(gb_sp_max, l, u, add = TRUE,
      col = 'purple', lwd = 2, lty = 6)
title(paste0('Sensitivity analysis + Spurious signal background:\n\n',
             'DATA: Cat',cat,'_',scenario,
             '_mu',mu_part_1,'p',mu_part_2, '; BKG: BW X Power0'))
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