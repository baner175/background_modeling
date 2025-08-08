rm(list = ls())

library(latex2exp)

# WITH BKG ONLY SAMPLE ##########################################

# BINNED; beta estimated:
B <- 1e4; n_phys <- 5e3; k <- 1e2
bkg_to_phys_ratio <- 2; eta_true <- 0

file_name <- paste0('Results/binned_test_eta_w_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')

# binned; beta known:
beta0 <- 3.87; B <- 1e4; n_phys <- 5e3; k <- 1e2
bkg_to_phys_ratio <- 2; eta_true <- 0

file_name <- paste0('Results/binned_test_eta_w_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')


# WITHOUT BKG ONLY SAMPLE ################################

# binned; beta estimated:
B <- 1e4; T_phys <- 5e3; k <- 1e2
lambda0 <- 0; eta_true <- 0

file_name <- paste0('Results/binned_test_eta_wo_bkg__',
                    'beta_estimated_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
                    'eta(',eta_true,')','.csv')

# binned; beta known:
beta0 <- 3.87; B <- 1e4; T_phys <- 5e3; k <- 1e2
lambda0 <- 0; eta_true <- 0

file_name <- paste0('Results/binned_test_eta_wo_bkg__',
                    'beta_known(',beta0,')_',
                    'B(',B,')_',
                    'T_phys(',T_phys,')_', 'n_bins(',k,')_',
                    'bkg_to_phys(',bkg_to_phys_ratio,')_',
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
