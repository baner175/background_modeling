library(latex2exp)

source('designing gb-2 bumps.R')

fb <- function(x) dtrunc(x, spec = 'norm',
                         mean = mean_back, sd = sd_back,
                         a = l, b = u)
f_mix <- function(x) eta_true*fs(x)+(1-eta_true)*fb(x)

xs <- seq(l, u, 0.01)
ys <- sapply(xs, f_mix)

My_Theme = theme(
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.text=element_text(size= 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot() + 
  geom_line(mapping = aes(x = xs, y = ys), linetype = 1,
                                 color = 'blue', linewidth = 1) + 
  ylab('f') + xlab('') + 
  theme_bw() + My_Theme +
  annotate('rect', fill = 'red', alpha = 0.4,
                                 xmax = M_upper, xmin = M_lower,
                                 ymax = Inf, ymin = -Inf) +
  annotate('text', x = c(2,3.5, 4.8),
           y = c(0.4,0.4,0.4), 
           label = c('Control Region (C)',
                     'Signal Region (S)',
                     'Control Region (C)'),
           size = 4)
