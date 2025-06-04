rm(list = ls())
library(dplyr)

### LOADING DATA ###############################################################

### change category, scenario and mu here ######################################
cat <- 0
scenario <- 'HLHC'
mu <- 1

mu_part_1 <- floor(mu); mu_part_2 <- (10*mu)%%10

file_name <- paste0('benchmark_toys/Cat',cat,'_',
                    scenario,'_mu',mu_part_1,'p',mu_part_2,'.csv')

obs <- read.csv(file_name, header = FALSE)[,1]

truth_file_name <- paste0('toymodels/h_Cat',cat,'_',
                          scenario,'_mu',mu_part_1,'.',mu_part_2,
                          '_mass__m.csv')

truth <- read.csv(truth_file_name, header = TRUE)

u <- 160; l <- 110
# dividing the range into 500 bins
bins <- seq(l, u, length.out = 501)
xi_truth <- (bins[1:500] + bins[2:501])/2
truth$xi <- xi_truth

# binning the data
k <- 100
bins <- seq(l, u, length.out = k+1)
bw <- (u-l)/k
ni <- sapply(1:k, function(i){
  sum((obs>bins[i])&(obs<=bins[i+1]))
})
N <- sum(ni)
xi <- (bins[1:k] + bins[2:(k+1)])/2
intervals <- lapply(xi, function(x) c(x + c(-1,1)*bw/2))
xi_collapsed <- sapply(xi_truth,
                       function(x) sapply(intervals,
                                          function(w) (x>=w[1]&x<w[2]))) %>% 
  apply(2, function(bool) xi[bool])
truth$xi_collapsed <- as.factor(xi_collapsed)


true_counts_df <- truth %>% group_by(xi_collapsed) %>% 
  summarise(collapsed_counts = sum(count_div10)) %>% as.data.frame()


hist(obs, breaks = bins, freq = TRUE, main = 'Comparing truth with realized data')
lines(x = xi, y = true_counts_df[,2], col = 'blue', lwd = 2)
