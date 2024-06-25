library(truncdist)
library(VGAM)
library(ggplot2)
library(latex2exp)

l <- 1; u <- 5
rate_back <- 3
rate_sig <- 1 #or 
shape_sig <- 0.01
rate_gb <- 2
eta <- 0.3

fb <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_back)
fs <- function(t) dtrunc(t, spec = 'pareto', a = l, b = u,
                         shape = shape_sig, scale = l)
gb <- function(t) dtrunc(t, spec = 'pareto', a = l, b = u,
                         shape = 50*shape_sig, scale = l)
# gb <- function(t) dtrunc(t, spec = 'exp', a = l, b = u, rate = rate_gb)

f_mix <- function(t) (1-eta)*fb(t) + eta*fs(t)

# delta:
integrate(function(t) (fs(t)/gb(t) - 1)*fb(t), l, u)


r = uniroot(function(t) fs(t) - f_mix(t), c(l, u))$root


My_Theme = theme(
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  legend.text=element_text(size= 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

ggplot(data.frame(x = c(l,u)), aes(x)) + 
  stat_function(fun = fb, 
                aes(linetype = 'true background', 
                    color = 'true background'), linewidth = 1) +
  stat_function(fun = fs,
                aes(linetype = 'signal',
                    color = 'signal'), linewidth = 1) +
  stat_function(fun = f_mix,
                aes(linetype = 'true mixture',
                    color = 'true mixture'), linewidth = 1) + 
  stat_function(fun = gb,
                aes(linetype = 'proposed background',
                    color = 'proposed background'), linewidth = 1) + 
  scale_linetype_manual('', values = c('true background' = 2,
                                       'true mixture' = 1,
                                       'proposed background' = 6,
                                       'signal' = 4))+
  scale_color_manual('',values = c('true background' = 'red',
                                'true mixture' = 'blue',
                                'proposed background' = 'orange',
                                'signal' = 'cyan')) +
  guides(linetype = guide_legend(title = ""),
         color = guide_legend(title = ""))+
  geom_vline(xintercept = r, alpha = 0.3,
             lwd = 1) +
  ylab('') + xlab('') + 
  annotate('text', x = c(r-0.5, 3),
           y = c(0.2,2),
           label = c(TeX('$\\R_{f_s < g_b}$'), TeX('$R_{f_s > g_b}$')),
           size = 8) +
  theme_bw() + My_Theme
