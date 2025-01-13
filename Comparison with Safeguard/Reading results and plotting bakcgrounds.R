# Plotting all the candidate backgrounds:

curve(gb_sf, lwd = 2.2, col = 'cyan', lty = 1, from = 1.15, to = 1.45,
      ylab = '', xlab = '')
curve(gb_test(x, fs_prop = 0), col = 'red',
      lty = 2, lwd = 2.2, add = TRUE)
curve(gb_test(x, fs_prop = 0.002), col = 'orange',
      lty = 3, lwd = 2.2, add = TRUE)
curve(gb_test(x, fs_prop = 0.005), col = 'brown',
      lty = 4, lwd = 2.2, add = TRUE)
curve(gb_test(x, fs_prop = 0.007), col = 'blue',
      lty = 5, lwd = 2.2, add = TRUE)
curve(gb_test(x, fs_prop = 0.01), col = 'green',
      lty = 6, lwd = 2.2, add = TRUE)
curve(fb_true, col = 'black',
      lty = 1, lwd = 2.2, add = TRUE)
curve(fs, col = 'skyblue',
      lty = 6, lwd = 2.2, add = TRUE)
abline(v = c(M_lower, M_upper), lwd = 2.2)

lambda_seq <- c(0, 0.002, 0.005, 0.007, 0.01)
legend(x = 1.35, y = 2.2,
       legend = c('safeguard',
                  latex2exp::TeX(sprintf(r'($g_b (\lambda = %f)$)', lambda_seq))),
       col = c('cyan', 'red', 'orange', 'brown', 'blue', 'green'),
       bty = 'n',
       lwd = 2.2,
       lty = 1:6,
       cex = 1.2)

# Reading files for type I error:
p_vals <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) p-values.csv', header = TRUE)
eta_hat <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) eta_hat.csv', header = TRUE)
test_stat <- read.csv('Simulation Results/Type 1 error (Ns+Nb=N) test_stat.csv', header = TRUE)

# Reading files for power:
eta_true <- 0.03
p_vals <- read.csv(paste0('Simulation Results/Power eta_', eta_true, ' (Ns+Nb=N) p-values.csv'), header = TRUE)
eta_hat <- read.csv(paste0('Simulation Results/Power eta_', eta_true, ' (Ns+Nb=N) eta_hat.csv'), header = TRUE)
test_stat <- read.csv(paste0('Simulation Results/Power eta_', eta_true, ' (Ns+Nb=N) test_stat.csv'), header = TRUE)


# Safeguard: 
mean(p_vals$Safeguard<=0.05) # corresponding delta: 0.01340076

# gb with fs_prop = 0.002:
mean(p_vals$fs_prop.0.002<=0.05) # corresponding delta: 0.016249060

# gb with fs_prop = 0.005:
mean(p_vals$fs_prop.0.005<=0.05) # corresponding delta: 0.007890316

# gb with fs_prop = 0.007:
mean(p_vals$fs_prop.0.007<=0.05) # corresponding delta: 0.002424086

# gb with fs_prop = 0.01:
mean(p_vals$fs_prop.0.01<=0.05) # corresponding delta: -0.005624272


# Reading results for binned data:

# Reading files for power:
eta_true <- 0
p_vals <- read.csv(paste0('Simulation Results/Binned - Power eta_', eta_true, ' (Ns+Nb=N) p-values.csv'), header = TRUE)
eta_hat <- read.csv(paste0('Simulation Results/Binned - Power eta_', eta_true, ' (Ns+Nb=N) eta_hat.csv'), header = TRUE)
test_stat <- read.csv(paste0('Simulation Results/Binned - Power eta_', eta_true, ' (Ns+Nb=N) test_stat.csv'), header = TRUE)

# gb with fs_prop = 0.002:
mean(p_vals$fs_prop.0.002<=0.05) # corresponding delta: 0.016249060

# gb with fs_prop = 0.005:
mean(p_vals$fs_prop.0.005<=0.05) # corresponding delta: 0.007890316

# gb with fs_prop = 0.007:
mean(p_vals$fs_prop.0.007<=0.05) # corresponding delta: 0.002424086

# gb with fs_prop = 0.01:
mean(p_vals$fs_prop.0.01<=0.05) # corresponding delta: -0.005624272

