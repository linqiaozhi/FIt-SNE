# Note: this script should be run after setting the FIt-SNE directory as the current one (e.g. setwd())

source('fast_tsne.R')

require(MASS);
N <- 1E4;
d <- 3;
input_data <- rbind(mvrnorm(n = N/2, rep(0, d), diag(d)),
	                  mvrnorm(n = N/2, rep(100, d), diag(d)))
Y2 <- fftRtsne(input_data, 2, max_iter = 50 );
plot(Y2[,1],Y2[,2],col=c(rep(1,N/2), rep(2,N/2))) # Plot the result

