source('fast_tsne.R', chdir=T)

library(gridExtra)
require(MASS);
library(FNN)
library(ggplot2)



d <- 3
N <- 200
input_data <- rbind(matrix( rep(1,N/2 * d), ncol=d), matrix(rep(2,N/2 * d), ncol=d))
dfs <- c(-0.4, -0.2, -0.1, 0.1,0.25, 0.5, 0.75, 1)
dfs <- c(0.001,0.01, 0.1,0.25, 0.5, 0.75, 1)
dfs <- c( .1)
col <- c(rep(1,N/2), rep(2,N/2))
p <- list()
max_iter = 600
for(i in 1:length(dfs)){
    Yg <- fftRtsne(input_data, 2,  df=dfs[i] ,ann_not_vptree=FALSE, max_iter =max_iter,rand_seed=1, sigma=2,  perplexity=-1, theta=0, learning_rate=1, exaggeration_factor=1);
    Yg_dist <-  min(get.knnx(Yg[1:(N/2),],Yg[( N/2 +1): N,], k=1)$nn.dist)
    p[[i]] <- ggplot(data.frame(Yg,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(show.legend = FALSE)  + ggtitle(sprintf('%d iters, df: %.4f, dist: %.2f', max_iter,dfs[i], Yg_dist)) + theme(plot.title = element_text(size=8))
}
do.call(grid.arrange,c(p,nrow=2))
(po =  exp(-1*sqrt(3)^2/(2*2^2)))
ps =  1
(po1 = po/( N*(po*(N/2)+ ps*(N/2-1))) )
(ps1 = ps/( N*(po*(N/2)+ ps*(N/2-1))) )




#n <- N*(N-1)/2
#n*(po1+ps1)
#
#n*(5.974406e-06 + 6.558386e-06)
#(epsilon = exp(-1*sqrt(3)^2/(2^2)))
#(p_ij = 2*epsilon/(N^2))
#(epsilon  = 2.048810e-05/2.981000e-05)
#sqrt((2/(epsilon) - 1)^(2/(dfs +1)) -1)
#
#
### Two points
#
#d <- 3
#N <- 200
#input_data <- rbind(matrix( rep(1,N/2 * d), ncol=d), matrix(rep(2,N/2 * d), ncol=d))
#dfs <- c(-0.4, -0.2, -0.1, 0.1,0.25, 0.5, 0.75, 1)
#dfs <- c(0.001,0.01, 0.1,0.25, 0.5, 0.75, 1)
#col <- c(rep(1,N/2), rep(2,N/2))
#p <- list()
#max_iter = 600
#for(i in 1:length(dfs)){
#    Yg <- fftRtsne(input_data, 2,  df=dfs[i] ,ann_not_vptree=FALSE, max_iter =max_iter,rand_seed=4, sigma=2,  perplexity=-1, theta=0, learning_rate=1, exaggeration_factor=1);
#    Yg_dist <-  min(get.knnx(Yg[1:(N/2),],Yg[( N/2 +1): N,], k=1)$nn.dist)
#    p[[i]] <- ggplot(data.frame(Yg,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(show.legend = FALSE)  + ggtitle(sprintf('%d iters, df: %.4f, dist: %.2f', max_iter,dfs[i], Yg_dist)) + theme(plot.title = element_text(size=8))
#}
#do.call(grid.arrange,c(p,nrow=2))
#
#
#
### With Gaussian, and perplexity
#
#d <- 3
#N <- 200
#set.seed(3)
#input_data <- rbind(mvrnorm(n = N/2, rep(0, d), diag(d)),mvrnorm(n = N/2, rep(100, d), diag(d)) )
#dfs <- c(0.01, 0.1,0.25, 0.5, 0.75, 1)
#dfs <- c( 1)
##dfs <- c(-0.2, -0.1, 0.001,0.01, 0.1,0.25, 0.5, 0.75, 1)
#col <- c(rep(1,N/2),rep(2,N/2))
#p <- list()
#max_iter = 600
#for(i in 1:length(dfs)){
#    #Yg <- fftRtsne(input_data, 2,  df=dfs[i] ,ann_not_vptree=FALSE, max_iter =max_iter,rand_seed=4, sigma=2,  perplexity=-1, theta=0, learning_rate=1, exaggeration_factor=1);
#    Yg <- fftRtsne(input_data, 2,  df=dfs[i] ,ann_not_vptree=FALSE, max_iter =max_iter,rand_seed=3, learning_rate=1, theta=0);
#    Yg_dist <-  min(get.knnx(Yg[1:(N/2),],Yg[( N/2 +1): N,], k=1)$nn.dist)
#    p[[i]] <- ggplot(data.frame(Yg,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(show.legend = FALSE)  + ggtitle(sprintf('%d iters, df: %.4f, dist: %.2f', max_iter,dfs[i], Yg_dist)) + theme(plot.title = element_text(size=8))
#}
#do.call(grid.arrange,c(p,nrow=2))
#
#
#p <- list()
#for(i in 1:length(dfs)){
#    Yg <- fftRtsne(input_data, 2,  df=dfs[i] ,ann_not_vptree=FALSE, max_iter =max_iter,rand_seed=4,   theta=0, learning_rate=1, exaggeration_factor=1);
#    Yg_dist <-  min(get.knnx(Yg[1:(N/2),],Yg[( N/2 +1): N,], k=1)$nn.dist)
#    p[[i]] <- ggplot(data.frame(Yg,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(show.legend = FALSE)  + ggtitle(sprintf('%d iters, df: %.4f, dist: %.2f', max_iter,dfs[i], Yg_dist)) + theme(plot.title = element_text(size=8))
#}
#do.call(grid.arrange,c(p,nrow=2))
#
#
#
#
#
#
#
#Ns <- c(1E3, 5E3, 1E4,2.5E4, 5E4);
#Ys <- list();
#
#N <- 90;
#
#d <- 10;
#input_data <- rbind(matrix( rep(1,N/2 * d), ncol=d), matrix(rep(2,N/2 * d), ncol=d))
#col <- c(rep(1,N/2), rep(2,N/2))
#Y <- fftRtsne(input_data, ann_not_vptree=FALSE, max_iter =500,rand_seed=3, sigma=1, perplexity=-1);
#
#(g0 <- ggplot(data.frame(Y,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() )
#
#
#
#
#Y2 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#Y3 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =2000,rand_seed=3);
#Y4 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =3000,rand_seed=3);
#Y5 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =4000,rand_seed=3);
#Y6 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =5000,rand_seed=3);
#
#Y2_dist <- nns <- min(get.knnx(Y2[1:(N/2),],Y2[( N/2 +1): N,], k=1)$nn.dist)
#Y3_dist <- nns <- min(get.knnx(Y3[1:(N/2),],Y3[( N/2 +1): N,], k=1)$nn.dist)
#Y4_dist <- nns <- min(get.knnx(Y4[1:(N/2),],Y4[( N/2 +1): N,], k=1)$nn.dist)
#Y5_dist <- nns <- min(get.knnx(Y5[1:(N/2),],Y5[( N/2 +1): N,], k=1)$nn.dist)
#Y6_dist <- nns <- min(get.knnx(Y6[1:(N/2),],Y6[( N/2 +1): N,], k=1)$nn.dist)
#g0 <- ggplot(data.frame(Y2,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle(sprintf('1000 iters: dist: %.2f', Y2_dist))
#g1 <- ggplot(data.frame(Y3,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle(sprintf('2000 iters: dist: %.2f', Y3_dist))
#g2 <- ggplot(data.frame(Y4,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle(sprintf('3000 iters: dist: %.2f', Y4_dist))
#g3 <- ggplot(data.frame(Y5,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle(sprintf('4000 iters: dist: %.2f', Y5_dist))
#g4 <- ggplot(data.frame(Y6,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle(sprintf('5000 iters: dist: %.2f', Y6_dist))
#gridExtra::grid.arrange(g0,g1,g2,g3,g4,nrow=1)
#
#df <- data.frame(Y3,col=factor(col))
#g1 <- ggplot(df, aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle('df=0.1')
##df <- data.frame(Y4,col=factor(col))
##g2 <- ggplot(df, aes(x=X1,y=X2,color=col)) + geom_point() 
#gridExtra::grid.arrange(g0,g1,nrow=1)
#
#
#
#
#
#
#
#
## Overshooting
#Y__ <- fftRtsne(input_data, 2,  df=0.01 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#g1 <- ggplot(data.frame(Y__,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=0.01, lr=200')
#Y__2 <- fftRtsne(input_data, 2,  df=0.01 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3, learning_rate =1);
#g2 <- ggplot(data.frame(Y__2,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=0.01, lr=1')
#gridExtra::grid.arrange(g1,g2,nrow=1)
#
#
#
#gauss_input_data <- rbind(mvrnorm(n = N, rep(0, d), diag(d)))
#Yg <- fftRtsne(gauss_input_data, 2,  df=1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#g1 <- ggplot(data.frame(Yg), aes(x=X1,y=X2)) + geom_point() + ggtitle('df=1, lr=200')
#Yg2 <- fftRtsne(gauss_input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#g2 <- ggplot(data.frame(Yg2), aes(x=X1,y=X2)) + geom_point() + ggtitle('df=0.2, lr=200')
#gridExtra::grid.arrange(g1,g2,nrow=1)
#
#
#
#d <- 3
#N <- 500
#gauss_input_data <- rbind(matrix( rep(1,N/2 * d), ncol=d), matrix(rep(2,N/2 * d), ncol=d))
#gauss_input_data <- rbind(mvrnorm(n = N/2, rep(0, d), diag(d)),mvrnorm(n = N/2, rep(2, d), diag(d)) )
#
#Yg <- fftRtsne(gauss_input_data, 2,  df=1 ,ann_not_vptree=FALSE, max_iter =1000,rand_seed=3, sigma=2,  perplexity=-1, theta=0);
#min(get.knnx(Yg[1:(N/2),],Yg[( N/2 +1): N,], k=1)$nn.dist)
#
#Yg <- fftRtsne(gauss_input_data, 2,  df=1 ,ann_not_vptree=FALSE, max_iter =1000,rand_seed=3, theta=0);
#col <- c(rep(1,N/2), rep(2,N/2))
#ggplot(data.frame(Yg, color=col), aes(x=X1,y=X2, color=color)) + geom_point() + ggtitle('df=1, lr=200')
#
#
#
#
### Blocks
#
#N <- 5E4;
#d <- 3;
#epsilon <- 1.08
#input_data <- rbind(matrix( runif(N*d, min=0,max=1), ncol=d))
#input_data[1:(N/2),1] <- input_data[1:(N/2),1]+ epsilon
#col <- c(rep(1,N/2), rep(2,N/2))
#Y <- fftRtsne(input_data, 2,  df=1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#Y2 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#g <- ggplot(data.frame(Y,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=1, lr=1')
#g2 <- ggplot(data.frame(Y2,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=1, lr=1')
#gridExtra::grid.arrange(g,g2,nrow=1)
#
#
#n <- 10
#seqs <- seq(0,1,by=1/n);
#M <- plot3D::mesh(seqs, seqs,seqs)
#input_data <- rbind( cbind(M$x,M$y,M$z), cbind(M$x+5,M$y,M$z))
#Y <- fftRtsne(input_data, 2,  df=1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#Y2 <- fftRtsne(input_data, 2,  df=0.1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#
#col <- c(rep(1,nrow(input_data)/2), rep(2,nrow(input_data)/2))
#g <- ggplot(data.frame(Y,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=1, lr=1')
#g2 <- ggplot(data.frame(Y2,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point() + ggtitle('df=1, lr=1')
#gridExtra::grid.arrange(g,g2,nrow=1)
#
#
#
### Seperate
#
#N <- 5E4;
#d <- 3;
#epsilon <- 1.08
#input_data <- rbind(matrix( runif(N*d, min=0,max=1), ncol=d))
#input_data[1:(N/2),1] <- input_data[1:(N/2),1]+ epsilon
#col <- c(rep(1,N/2), rep(2,N/2))
#Y <- fftRtsne(input_data, 2,  df=1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3);
#
#
#
#
#
#
#
#
#
##input_data <- matrix(0, nrow=N,ncol=d)
##diag(input_data) <- 1
##input_data <- rbind(mvrnorm(n = N, rep(0, d), diag(d)))
##input_data2 <- input_data/sqrt(rowSums(input_data^2))
##input_data3 <- rbind(input_data2[(N/2+1):N,] + 0.05, input_data2[1:(N/2),]-0.05)
##input_data3 <- input_data
##plot(input_data3[,1], input_data3[,2],col=c(rep(1,N/2), rep(2,N/2)))
#
#Y2 <- fftRtsne(input_data, 2,  df=1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3, learning_rate=50, exaggeration_factor=1);
#
#Y3 <- fftRtsne(input_data, 2,  df=-0.1 , ann_not_vptree=FALSE, max_iter =1000,rand_seed=3, learning_rate=50, exaggeration_factor=1);
#
#
#df <- data.frame(Y2,col=factor(col))
#g0 <- ggplot(df, aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle('df=1')
#df <- data.frame(Y3,col=factor(col))
#g1 <- ggplot(df, aes(x=X1,y=X2,color=col)) + geom_point()  + ggtitle('df=0.1')
##df <- data.frame(Y4,col=factor(col))
##g2 <- ggplot(df, aes(x=X1,y=X2,color=col)) + geom_point() 
#gridExtra::grid.arrange(g0,g1,nrow=1)
#
#
#
#plot(Y2[,1],Y2[,2],col=c(rep(1,N/2), rep(2,N/2))) # Plot the result
#plot(Y3[,1],Y3[,2],col=c(rep(1,N/2), rep(2,N/2))) # Plot the result
#
#Qs2 <- read.csv('temp/Qs.csv', header =F);
#Qi2 <- Qs2[,1]
#Qs2 <- Qs2[,-1]
#Qs2[1,]
#
#Y3 <- fftRtsne(input_data3, 2,  df=0.1 , ann_not_vptree=FALSE);
#Qs2 <- read.csv('temp/Qs.csv', header =F);
#Qi2 <- Qs2[,1]
#Qs2 <- Qs2[,-1]
#max(Qs2[1,])
#
#k <- 20
#p <- 2
#nns <- get.knn(input_data3,k=k)
#col <- (rep(0.1, N))
#col[nns$nn.index[p,]] <-1:k +1
#col[p] <- 1
#
#df <- data.frame(Y2,col=factor(col),siz=col)
#g1 <- ggplot(df, aes(x=X1,y=X2,color=col,size=siz)) + geom_point() + scale_size(range=c(2,6))
#df <- data.frame(Y3,col=factor(col),siz=col)
#g2 <- ggplot(df, aes(x=X1,y=X2,color=col,size=siz)) + geom_point() + scale_size(range=c(2,6))
#gridExtra::grid.arrange(g1,g2,nrow=1)
#
#plot(Y2[,1],Y2[,2],col=c(rep(1,N/2), rep(2,N/2))) # Plot the result
#plot(Y3[,1],Y3[,2],col=c(rep(1,N/2), rep(2,N/2))) # Plot the result
