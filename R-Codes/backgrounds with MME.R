source('bkg_MLE functions.R')
source('functions.R')

library(ggplot2)

mean_sig <-  0.4
mean_back <- 5
sd_sig <- 0.1
sd_back <- 5

data <- make_data(mean_sig = 0.4, 
                  sd_sig = 0.1,
                  mean_back = 5,
                  sd_back = 5,
                  bkg_prop = 0.97,
                  seed = 12345)

data_back <- data$background
data <- data$observed

old_basis <- c(T1_s, T2_s, T3_s)
new_basis <- c(T1, T2, T3)

bkg_estimation_MME <- function(data, nfun, new = FALSE, ...)
{
  tau <- NULL
  T_mat <- NULL
  if(new){
    basis <- new_basis
  }else{
    basis <- old_basis
  }
  for (j in 1:nfun) {
    Tj <- basis[j][[1]]
    vec <- sapply(data, function(t)
    {
      Tj(t,...)
    })
    tau <- c(tau, mean(vec))
    T_mat <- rbind(T_mat, vec)
  }
  return(list(T_mat = T_mat, tau = tau))
}

xs <- seq(0,1,0.01)
y_real <- sapply(xs, actual, mean_sig = 0.4, 
                 sd_sig = 0.1,
                 mean_back = 5,
                 sd_back = 5,
                 bkg_prop = 0.97)
# one function

res_mixD_oldB_1 <- bkg_estimation_MME(data = data,
                                    nfun = 1,
                                    new = FALSE)

den_mixD_oldB_1 <- apply(res_mixD_oldB_1$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_oldB_1$tau))
})


res_backD_oldB_1 <- bkg_estimation_MME(data = data_back,
                                    nfun = 1,
                                    new = FALSE)
den_backD_oldB_1 <- apply(res_backD_oldB_1$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_oldB_1$tau))
})

res_mixD_newB_1 <- bkg_estimation_MME(data = data,
                                    nfun = 1,
                                    new = TRUE,
                                    mean = mean_sig,
                                    sd = sd_sig)

den_mixD_newB_1 <- apply(res_mixD_newB_1$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_newB_1$tau))
})

res_backD_newB_1 <- bkg_estimation_MME(data = data_back,
                                    nfun = 1,
                                    new = TRUE,
                                    mean = mean_sig,
                                    sd = sd_sig)
den_backD_newB_1 <- apply(res_backD_newB_1$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_newB_1$tau))
})


plt_T1 <- ggplot(mapping=aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = data, y = den_mixD_newB_1,
                          color = 'mixed data, new basis'),
            linetype = 'dashed', lwd = 1.1, alpha = 0.3) +
  geom_line(mapping = aes(x = data_back, y = den_backD_newB_1,
                          color = 'background data, new basis'),
            linetype = 'dashed', lwd = 1.1,  alpha = 0.3) + 
  geom_line(mapping = aes(x = data_back, y = den_backD_oldB_1,
                          color = 'background data, old basis'),
            linetype = 'dashed', lwd = 1.1,  alpha = 0.3) +
  geom_line(mapping = aes(x = data, y = den_mixD_oldB_1,
                          color = 'mixed data, old basis'),
            linetype = 'dashed', lwd = 1.1,  alpha = 0.3) + 
  geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                 a = 0, b = 1),
                          color = 'actual background')
            , lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dunif(xs),
                          color = 'proposed background'), lwd = 1.1, 
            linetype = 'dashed', alpha = 0.3) +
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) + 
  labs(x = "background data",
       y = "",
       color = "Legend") +
  scale_color_manual(values = c('mixed data, new basis' = 'green',
                                'background data, new basis' = 'orange',
                                'background data, old basis' = 'blue',
                                'mixed data, old basis' = 'purple',
                                'actual background' = 'black',
                                'proposed background' = 'red',
                                'true model' = 'pink')) + 
  ggtitle(bquote(T[1]))

#------------------------------------------------------------------------

res_mixD_oldB_2 <- bkg_estimation_MME(data = data,
                                      nfun = 2,
                                      new = FALSE)

den_mixD_oldB_2 <- apply(res_mixD_oldB_2$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_oldB_2$tau))
})


res_backD_oldB_2 <- bkg_estimation_MME(data = data_back,
                                       nfun = 2,
                                       new = FALSE)
den_backD_oldB_2 <- apply(res_backD_oldB_2$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_oldB_2$tau))
})

res_mixD_newB_2 <- bkg_estimation_MME(data = data,
                                      nfun = 2,
                                      new = TRUE,
                                      mean = mean_sig,
                                      sd = sd_sig)

den_mixD_newB_2 <- apply(res_mixD_newB_2$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_newB_2$tau))
})

res_backD_newB_2 <- bkg_estimation_MME(data = data_back,
                                       nfun = 2,
                                       new = TRUE,
                                       mean = mean_sig,
                                       sd = sd_sig)
den_backD_newB_2 <- apply(res_backD_newB_2$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_newB_2$tau))
})

plt_T2 <- ggplot(mapping=aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = data, y = den_mixD_newB_2,
                          color = 'mixed data, new basis'),
            linetype = 'dashed', lwd = 1.1, alpha = 0.3) +
  geom_line(mapping = aes(x = data_back, y = den_backD_newB_2,
                          color = 'background data, new basis'),
            linetype = 'dashed', lwd = 1.1, alpha = 0.3) + 
  geom_line(mapping = aes(x = data_back, y = den_backD_oldB_2,
                          color = 'background data, old basis'),
            linetype = 'dashed', lwd = 1.1,alpha = 0.3) +
  geom_line(mapping = aes(x = data, y = den_mixD_oldB_2,
                          color = 'mixed data, old basis'),
            linetype = 'dashed', lwd = 1.1,alpha = 0.3) + 
  geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                 a = 0, b = 1),
                          color = 'actual background')
            , lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dunif(xs),
                          color = 'proposed background'), lwd = 1.1,
            alpha = 0.3) +
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) + 
  labs(x = "background data",
       y = "",
       color = "Legend") +
  scale_color_manual(values = c('mixed data, new basis' = 'green',
                                'background data, new basis' = 'orange',
                                'background data, old basis' = 'blue',
                                'mixed data, old basis' = 'purple',
                                'actual background' = 'black',
                                'proposed background' = 'red',
                                'true model' = 'pink')) + 
  ggtitle(bquote(T[1] + T[2]))
#------------------------------------------------------------------------

res_mixD_oldB_3 <- bkg_estimation_MME(data = data,
                                      nfun = 3,
                                      new = FALSE)

den_mixD_oldB_3 <- apply(res_mixD_oldB_3$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_oldB_3$tau))
})


res_backD_oldB_3 <- bkg_estimation_MME(data = data_back,
                                       nfun = 3,
                                       new = FALSE)
den_backD_oldB_3 <- apply(res_backD_oldB_3$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_oldB_3$tau))
})

res_mixD_newB_3 <- bkg_estimation_MME(data = data,
                                      nfun = 3,
                                      new = TRUE,
                                      mean = mean_sig,
                                      sd = sd_sig)

den_mixD_newB_3 <- apply(res_mixD_newB_3$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_mixD_newB_3$tau))
})

res_backD_newB_3 <- bkg_estimation_MME(data = data_back,
                                       nfun = 3,
                                       new = TRUE,
                                       mean = mean_sig,
                                       sd = sd_sig)
den_backD_newB_3 <- apply(res_backD_newB_3$T_mat, 2, function(t){
  1 + as.numeric(crossprod(t, res_backD_newB_3$tau))
})

plt_T3 <- ggplot(mapping=aes(x = data_back)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 bins = 20, fill = 'steelblue', col = 'black') + 
  geom_line(mapping = aes(x = data, y = den_mixD_newB_3,
                          color = 'mixed data, new basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = data_back, y = den_backD_newB_3,
                          color = 'background data, new basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = data_back, y = den_backD_oldB_3,
                          color = 'background data, old basis'),
            linetype = 'dashed', lwd = 1.1) +
  geom_line(mapping = aes(x = data, y = den_mixD_oldB_3,
                          color = 'mixed data, old basis'),
            linetype = 'dashed', lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dtruncnorm(xs, mean = 5, sd = 5,
                                                 a = 0, b = 1),
                          color = 'actual background')
            , lwd = 1.1) + 
  geom_line(mapping = aes(x = xs, y = dunif(xs),
                          color = 'proposed background'), lwd = 1.1,
            linetype = 'dashed') +
  geom_line(mapping = aes(x = xs, y = y_real,
                          color = 'true model'), lwd = 1.1) + 
  labs(x = "background data",
       y = "",
       color = "Legend") +
  scale_color_manual(values = c('mixed data, new basis' = 'green',
                                'background data, new basis' = 'orange',
                                'background data, old basis' = 'blue',
                                'mixed data, old basis' = 'purple',
                                'actual background' = 'black',
                                'proposed background' = 'red',
                                'true model' = 'pink')) + 
  ggtitle(bquote(T[1] + T[2] + T[3]))

ggpubr::ggarrange(plt_T1, plt_T2, plt_T3,
                  ncol = 3, nrow = 1,
                  common.legend = TRUE, legend = 'bottom')

