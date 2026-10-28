rm(list = ls())
library(truncdist)
library(VGAM)
library(doSNOW)
library(parallel)
library(foreach)
library(latex2exp)

npLRT_neg_log_lik <- function(eta, ni, mi, si){
  bi <- mi/sum(mi)
  return(-sum(ni*log(eta*si + (1-eta)*bi)))
}


simulated_power_binned_npLRT <- function(eta, nbins, T_phys,
                                         r, nsims, seed = 12345,
                                         signif.level = 0.05){
  cat("\n--------------------------------------\n")
  cat(sprintf('Evaluating empirical power at eta = %.2f', eta))
  cat("\n--------------------------------------\n")
  
  k <- nbins; B <- nsims; T_bkg <- r*T_phys
  bin_ends <- seq(l, u, length.out = k+1)
  xi <- (bin_ends[-1] + bin_ends[-(k+1)])/2
  si <- sapply(1:k, function(i){
    integrate(fs, bin_ends[i], bin_ends[i+1])$value
  })
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  clusterExport(cl, c('l', 'u', 'mean_sig', 'sd_sig', 'bkg_rate', 'bkg_shape',
                      'npLRT_neg_log_lik'))
  pb <- txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  test_stat_LRT <- foreach(i = 1:B, .combine = c,
                           .packages = c('truncdist', 'VGAM'),
                           .options.snow = opts) %dopar%
    {
      set.seed(seeds[i])
      # total sample sizes:
      M <- rpois(1, lambda = T_bkg); N <- rpois(1, lambda = T_phys)
      
      # bkg-only sample:
      bkg_samp <- rtrunc(M, a = l, b = u, spec = 'gamma',
                         rate = bkg_rate, shape = bkg_shape)
      
      # physics-sample:
      s_samp <- rtrunc(N, a = l, b = u, spec = 'norm',
                       mean = mean_sig, sd = sd_sig)
      b_samp <- rtrunc(N, a = l, b = u, spec = 'gamma',
                       rate = bkg_rate, shape = bkg_shape)
      u_mask <- runif(N)
      phys_samp <- ifelse(u_mask <= eta, s_samp, b_samp)
      
      ni <- sapply(1:k, function(i){
        sum((phys_samp>bin_ends[i])&(phys_samp<=bin_ends[i+1]))
      })
      mi <- sapply(1:k, function(i){
        sum((bkg_samp>bin_ends[i])&(bkg_samp<=bin_ends[i+1]))
      })
      mi_cens <- mi*(mi>=1) + (mi<1)
      opt <- nlminb(start = 0.01,
                    objective = npLRT_neg_log_lik,
                    lower = 0, upper = 1,
                    ni = ni, mi = mi_cens, si = si)
      eta_hat <- opt$par
      bi_hat <- mi_cens/sum(mi_cens)
      pi_eta_hat <- eta_hat*si + (1-eta_hat)*bi_hat
      
      2*sum(ni*(log(pi_eta_hat) - log(bi_hat)))
    }
  
  p_vals <- 0.5* pchisq(test_stat_LRT, df = 1, lower.tail = FALSE)
  simulated_power <- mean(p_vals < signif.level)
  return(list(test_stat = test_stat_LRT, power = simulated_power))
}

#parameters for the signal
mean_sig <- 1.28
sd_sig <- 0.02
l <- 1; u <- 2

# signal density
fs <- function(x, mean = mean_sig) dtrunc(x, a = l, b = u, spec = 'norm', mean = mean, sd = sd_sig)

#parameter for the true background
bkg_rate <- 3.3; bkg_shape <- 0.5

# true bkg density
fb_true <- function(x) dtrunc(x, a = l, b = u, spec = 'gamma',
                              rate = bkg_rate, shape = bkg_shape)


B <- 1e5; n_phys <- 2e3
bkg_to_phys_ratio <- 2; eta <- 0
n_bkg <- n_phys*bkg_to_phys_ratio


res_1e2 <- simulated_power_binned_npLRT(eta = eta, nbins = 1e2, T_phys = n_phys,
                                    r = bkg_to_phys_ratio, nsims = B,
                                    seed = 12345, signif.level = 0.05)
ecdf_LRT_1e2 <- ecdf(res_1e2$test_stat)

res_50 <- simulated_power_binned_npLRT(eta = eta, nbins = 50, T_phys = n_phys,
                                        r = bkg_to_phys_ratio, nsims = B,
                                        seed = 12345, signif.level = 0.05)
ecdf_LRT_50 <- ecdf(res_50$test_stat)

res_30 <- simulated_power_binned_npLRT(eta = eta, nbins = 30, T_phys = n_phys,
                                       r = bkg_to_phys_ratio, nsims = B,
                                       seed = 12345, signif.level = 0.05)
ecdf_LRT_30 <- ecdf(res_30$test_stat)

# takes about 10 minutes

op <- par(no.readonly = TRUE)
par(mgp = c(2, 0.8, 0))
curve(0.5+0.5*pchisq(x, df = 1), lwd = 4, col = 'black',
      lty = 1, from = 0, to = 20, ylim = c(0,1), ylab = TeX('$P(X \\leq x)$'),
      cex.lab = 2, cex.axis = 2)
abline(h = 0.5, v = 0, col = 'grey', lty = 2, lwd = 4)
curve(ecdf_LRT_1e2, ylim = c(0,1), lwd = 4, add = TRUE, col = 'orange', lty = 2)
curve(ecdf_LRT_50, ylim = c(0,1), lwd = 4, add = TRUE, col = 'red', lty = 6)
curve(ecdf_LRT_30, ylim = c(0,1), lwd = 4, add = TRUE, col = 'blue', lty = 4)

legend(x = 11.5, y = 0.55,
       legend = c(
         TeX('$\\bar{\\chi}^2_{0 1}$'),
         paste0(c(30,50,100), ' bins')
       ),
       col = c('black', 'blue', 'red', 'orange'),
       lty = c(1, 4, 6, 2),
       bty = 'n',
       lwd = 4,
       seg.len = 2,
       y.intersp = 1.5, cex = 2)
par(op)
