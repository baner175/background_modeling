# library(devtools)
# devtools::install_github(repo = "rohitpatra/mixmodel",ref='main')
library(mixmodel)

obs <- read.table('Data_ex1.txt')$x

est.default <- mix.model(obs, method = "fixed", gridsize = 600)

print(est.default) #"Estimate of alp is 0.995"

est.fixed <- mix.model(obs, method = "fixed",
                       c.n = .05*log(log(length(obs))),
                       gridsize = 600)
print(est.fixed) #"Estimate of alp is 0.996"
