rm(list = ls())
generate_unbinned <- function(file_path, seed = 12345,
                              lower, upper, bw, mu_s)
{
  l = lower; u = upper; bw = 0.1
  binned <- read.csv(file_path, header = TRUE)
  bins <- seq(l, u, bw)
  obs <- c()
  
  set.seed(seed = seed)
  for(i in 1:(length(bins))-1)
  {
    obs <- c(obs, runif(binned$Content[i], bins[i], bins[i+1]))
  }
  
  hist(obs, probability = TRUE, breaks = 50)
  abline(v = mu_s, col = 'red', lwd = 1.2)
  # kde <- kdensity::kdensity(obs, bw = bw)
  # y_kde <- sapply(obs, kde)
  # lines(obs, y_kde, col = 'black')
  
  return(obs)
}


# Toy Data A:

obs <- generate_unbinned(file_path = 'Toy_Data_A/Data/histogram_data_cat0.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_A/Data/pseudo_unbinned_cat0.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_A/Data/histogram_data_cat1.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_A/Data/pseudo_unbinned_cat1.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_A/Data/histogram_data_cat2.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_A/Data/pseudo_unbinned_cat2.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_A/Data/histogram_data_cat3.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_A/Data/pseudo_unbinned_cat3.csv", row.names = FALSE)


# Toy Data B:

obs <- generate_unbinned(file_path = 'Toy_Data_B/Data/histogram_data_cat0.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_B/Data/pseudo_unbinned_cat0.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_B/Data/histogram_data_cat1.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_B/Data/pseudo_unbinned_cat1.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_B/Data/histogram_data_cat2.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_B/Data/pseudo_unbinned_cat2.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_B/Data/histogram_data_cat3.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_B/Data/pseudo_unbinned_cat3.csv", row.names = FALSE)


# Toy Data C:

obs <- generate_unbinned(file_path = 'Toy_Data_C/Data/histogram_data_cat0.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_C/Data/pseudo_unbinned_cat0.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_C/Data/histogram_data_cat1.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_C/Data/pseudo_unbinned_cat1.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_C/Data/histogram_data_cat2.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_C/Data/pseudo_unbinned_cat2.csv", row.names = FALSE)

obs <- generate_unbinned(file_path = 'Toy_Data_C/Data/histogram_data_cat3.csv',
                         seed = 12345, 
                         lower = 110, upper = 160, bw = 0.1, mu_s = 125)
write.csv(x = obs, file = "Toy_Data_C/Data/pseudo_unbinned_cat3.csv", row.names = FALSE)






