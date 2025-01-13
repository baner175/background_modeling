rm(list = ls())

l <- 110; u <- 160

obs1 <- read.csv('Data/toy1.csv', header = TRUE)$x
hist(obs1, probability = TRUE, breaks = 50, main = 'Toy 1')
abline(v=125, lwd = 2, col = 'red')
obs_normalised <- scale(obs1)
kde_ <- kdensity::kdensity(obs_normalised)
kde_unnormed<- function(t) kde_((t-mean(obs1))/sd(obs1))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
curve(kde, l, u, col = 'black', lwd = 2, add = TRUE)

obs2 <- read.csv('Data/toy2.csv', header = TRUE)$x
hist(obs2, probability = TRUE, breaks = 50, main = 'Toy 2')
abline(v=125, lwd = 2, col = 'red')
obs_normalised <- scale(obs2)
kde_ <- kdensity::kdensity(obs_normalised)
kde_unnormed<- function(t) kde_((t-mean(obs2))/sd(obs2))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
curve(kde, l, u, col = 'black', lwd = 2, add = TRUE)

obs3 <- read.csv('Data/toy3.csv', header = TRUE)$x
hist(obs3, probability = TRUE, breaks = 50, main = 'Toy 3')
abline(v=125, lwd = 2, col = 'red')
obs_normalised <- scale(obs3)
kde_ <- kdensity::kdensity(obs_normalised)
kde_unnormed<- function(t) kde_((t-mean(obs3))/sd(obs3))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
curve(kde, l, u, col = 'black', lwd = 2, add = TRUE)


obs4 <- read.csv('Data/toy4.csv', header = TRUE)$x
hist(obs4, probability = TRUE, breaks = 50, main = 'Toy 4')
abline(v=125, lwd = 2, col = 'red')
obs_normalised <- scale(obs4)
kde_ <- kdensity::kdensity(obs_normalised)
kde_unnormed<- function(t) kde_((t-mean(obs4))/sd(obs4))
kde_tot <- integrate(kde_unnormed, l, u)$value
kde <- function(x) kde_unnormed(x)/kde_tot
curve(kde, l, u, col = 'black', lwd = 2, add = TRUE)

