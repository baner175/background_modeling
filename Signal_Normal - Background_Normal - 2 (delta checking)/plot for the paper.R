library(latex2exp)

source('designing gb-2 bumps.R')

fb <- function(x) dtrunc(x, spec = 'norm',
                         mean = mean_back, sd = sd_back,
                         a = l, b = u)
f_mix <- function(x) eta_true*fs(x)+(1-eta_true)*fb(x)

My_Theme = theme(
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.text=element_text(size= 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(data.frame(x = c(l,u)), aes(x)) + 
  stat_function(fun = fb, aes(linetype = 'true background',
                              color = 'true background'), linewidth = 1) + 
  stat_function(fun = f_mix, aes(linetype = 'true mixture',
                                 color = 'true mixture'), linewidth = 1) + 
  stat_function(fun = gb, aes(linetype = 'proposed background',
                              color = 'proposed background'), linewidth = 1) + 
  scale_linetype_manual('', values = c('true background' = 2,
                                       'true mixture' = 1,
                                       'proposed background' = 6))+
  scale_color_manual('',values = c('true background' = 'red',
                                   'true mixture' = 'blue',
                                   'proposed background' = 'orange')) +
  guides(linetype = guide_legend(title = ""),
         color = guide_legend(title = ""))+
  geom_vline(xintercept = c(M_lower, M_upper), alpha = 0.3,
             lwd = 1) + 
  ylab('') + xlab('') + 
  annotate('text', x = c(M_lower-0.1, M_upper+0.1),
           y = c(0,0), 
           label = c(TeX('$\\mu_s-d_\\epsilon$'), TeX('$\\mu_s+d_\\epsilon$')),
           size = 10) + 
  theme_bw() + My_Theme
