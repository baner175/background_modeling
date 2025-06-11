rm(list = ls())

library(latex2exp)

# WITH BKG ONLY SAMPLE ##########################################

# UNBINNED; beta estimated:
B <- 1e4; n_phys <- 5e3
bkg_to_phys_ratio <- 2; eta_true <- 0

file_name <- paste0('Results/unbinned_test_eta_w_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')

# UNBINNED; beta known:
beta0 <- 3.87; B <- 1e4; n_phys <- 5e3
bkg_to_phys_ratio <- 2; eta_true <- 0

file_name <- paste0('Results/unbinned_test_eta_w_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'n_phys(',n_phys,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')', '.csv')

# WITHOUT BKG ONLY SAMPLE ################################

# UNBINNED; beta estimated:
B <- 1e4; n_samp <- 5e3
lambda0 <- 0; eta_true <- 0

file_name <- paste0('Results/unbinned_test_eta_wo_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'n_samp(',n_samp,')_',
                    'lambda(',lambda0,')_',
                    'eta(',eta_true,')','.csv')

# UNBINNED; beta known:
beta0 <- 3.87; B <- 1e4; n_samp <- 5e3
lambda0 <- 0; eta_true <- 0

file_name <- paste0('Results/unbinned_test_eta_wo_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'n_samp(',n_samp,')_',
                    'lambda(',lambda0,')_',
                    'eta(',eta_true,')','.csv')

# plotting ecdf for eta test statistic and comparing with pnorm
df <- read.csv(file_name, header = TRUE)
plt_title <- strsplit(file_name,'/')[[1]][2]
plt_title <- strsplit(plt_title, '.csv')[[1]]
plt_title <-  paste0(strsplit(plt_title, '__beta')[[1]], 
                     collapse = '\nbeta')
hist(df$test_stat, probability = TRUE, 
     breaks = 50, main = '')
curve(dnorm, add = TRUE, lwd = 2)
title(main = plt_title)

Fn <- ecdf(df$test_stat)
curve(Fn, from = -2, to = 3, col = 'red', lwd = 2, lty = 1)
curve(pnorm, add = TRUE, col = 'black', lwd = 2, lty = 2)
legend('bottomright', 
       legend = c('std normal CDF', 'Empirical CDF'),
       col = c('black', 'red'), lty = c(2,1), lwd = 2,
       bty = 'n')
title(main = plt_title)
print(mean(df$test_stat>qnorm(0.05, lower.tail = FALSE)))

####################################################################
# LRT simulation:
B <- 1e4; n_samp <- 1e3; eta_true <- 0
file_name <- paste0('Results/LRT', 
                    '_B_', B, 
                    '_n_samp_', n_samp,
                    '_eta_', eta_true,
                    '.csv')
df <- read.csv(file_name, header = TRUE)

if(eta_true == 0){
  Fn_star <- ecdf(df[,1])
  Fn_0 <- ecdf(df[,2])
  theo_CDF_0 <- function(x) (0.5 + 0.5*pchisq(x, df = 1))*(x>=0)
  theo_CDF_1 <- function(x) pchisq(x, df = 1)
  curve(theo_CDF_0, from = 0, to = 20, ylab = '',
        lwd = 2, lty = 1, col = 'black')
  curve(Fn_0, add = TRUE, lty = 2, col = 'red')
  curve(Fn_star, add = TRUE, lwd = 2, lty = 3, col = 'blue')
  curve(theo_CDF_1, add = TRUE, lwd = 2, lty = 3, col = 'brown')
  legend('bottomright', 
         legend = c(TeX('$\\frac{1}{2}\\delta_0 + \\frac{1}{2}\\chi^2_1$'),
                    TeX('Emp. CDF of $-2 log(\\lambda(0))$'),
                    TeX('Emp. CDF of $-2 log(\\lambda(\\tilde{\\eta_*}))$')),
         col = c('black', 'red', 'blue'), lty = 1:3,
         lwd = 2,
         bty = 'n')
  
}else{
  Fn_star <- ecdf(df[,1])
  theo_CDF <- function(x) pchisq(x, df = 1)
  curve(theo_CDF, from = 0, to = 20, ylab = '',
        lwd = 2, lty = 1, col = 'black')
  curve(Fn_star, add = TRUE, lty = 2, col = 'red')
  legend('bottomright', 
         legend = c(TeX('$\\chi^2_1$'),
                    TeX('Emp. CDF of $-2 log(\\lambda(\\tilde{\\eta_*}))$')),
         col = c('black', 'red'), lty = c(1,2),
         lwd = 2,
         bty = 'n')
}
