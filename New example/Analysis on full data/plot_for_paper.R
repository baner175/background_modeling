source('designing gb with 2 bumps.R')

mean_back <- 0.5; sd_back <- 2.5
eta_true <- 0.03
f_back <- function(t) dtrunc(t, spec = 'norm', a = l, b = u,
                             mean = mean_back, sd = sd_back)
f_mix <- function(t) eta_true*fs(t) + (1-eta_true)*f_back(t)

xs <- seq(l,u, 0.01)  
y_mix <- sapply(xs, f_mix)
y_back <- sapply(xs, f_back)
y_gb <- sapply(xs, gb)

txt_size = 20
My_Theme = theme(
  axis.title.x = element_text(size = txt_size),
  axis.text.x = element_text(size = txt_size),
  axis.title.y = element_text(size = txt_size),
  axis.text.y = element_text(size = txt_size),
  legend.text = element_text(size = txt_size),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(mapping = aes(x = obs)) + 
  geom_histogram(mapping = aes(y = after_stat(density)),
                 fill = 'white', col = 'black',
                 breaks = hs$breaks) + 
  geom_line(mapping = aes(x = xs, y = y_gb,
                          linetype = 'proposal background',
                          color = 'proposal background'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = y_back,
                          linetype = 'true background',
                          color = 'true background'),
            lwd = 1.2) + 
  geom_line(mapping = aes(x = xs, y = y_mix,
                          linetype = 'mixture density',
                          color = 'mixture density'),
            lwd = 1.2)  + 
  xlab('log(y)') + ylab('Density') +
  scale_color_manual('', 
                     values = c('proposal background' = 'orange',
                                'true background' = 'red',
                                'mixture density' = 'blue')) + 
  scale_linetype_manual('',
                        values = c('proposal background' = 6,
                                   'true background' = 2,
                                   'mixture density' = 1)) + 
  guides(linetype = guide_legend(title = ""),
         color = guide_legend(title = "")) +
  theme_bw() +
  My_Theme

