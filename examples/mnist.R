library(dbscan)
library(rsvd)
library(gridExtra)
library(ggplot2)
source('../fast_tsne.R', chdir=T)

show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1],yaxt='n',xaxt='n', ann=FALSE, col=col, ...)
}

load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('~/Downloads/mnist/train-images-idx3-ubyte')
  train$y <<- load_label_file('~/Downloads/mnist/train-labels-idx1-ubyte')
}

load_mnist()
X = train$x/255
y = train$y

# PCA
Xm <- colMeans(X)
Xc <- sweep(X, 2,Xm)
set.seed(3)
fastDecomp <- rsvd(Xc, 50,q=4);
fastDecomp$d
PCs<- fastDecomp$u %*% diag(fastDecomp$d);

Xg_all_0.1<- fftRtsne(PCs, 2,  df=0.55 ,rand_seed=4, learning_rate=100, max_iter=2000);
idx <- which(y == 0)
Xg_0_0.1 <- fftRtsne(PCs[idx,], 2,  df=0.55 ,rand_seed=4, learning_rate=100, max_iter=2000);
Xg_0_m0.9 <- fftRtsne(PCs[idx,], 2,  df=0.05 ,rand_seed=4, learning_rate=100, max_iter=2000);

(clust <- dbscan(Xg_all_0.1[idx,], eps= 2, minPts=100))
(clust2 <- dbscan(Xg_0_m0.9, eps= 5, minPts= 100))

t <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
col <- clust2$cluster

g1 <- ggplot(data.frame(Xg_all_0.1[idx,],col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(size=0.5)  +  theme(plot.title = element_text(size=8)) +  
    ggtitle(expression(paste(bold("A) "),"Embedding all, showing digit 0, df=0.1")))+ guides(color=FALSE) #+ t
g2 <- ggplot(data.frame(Xg_0_0.1,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(size=0.5)  +  theme(plot.title = element_text(size=8)) + 
    ggtitle(expression(paste(bold("B) "),"Embedding only 0, df=0.1")))+ guides(color=FALSE)# +t 
g3 <- ggplot(data.frame(Xg_0_m0.9,col=factor(col)), aes(x=X1,y=X2,color=col)) + geom_point(size=0.5)  +  theme(plot.title = element_text(size=8)) + 
    ggtitle(expression(paste(bold("C) "),"Embedding only 0, df=-0.9"))) + guides(color=FALSE)# +t
(g <- grid.arrange(g1,g2,g3, nrow=1))

ggsave("embeddings_single_digit.pdf", g)


par(mfrow=c(2,6),mar=c(1,1,1,1) )
for(i in 1:max(clust2$cluster)){
    meanify <- idx[sample(1:length(idx), 400)]
    meanify <- idx[clust2$cluster == i]
    clustmean <- colMeans(X[meanify,])
    show_digit(clustmean)
    title(sprintf("N=%d", length(meanify)))
}

dev.print(pdf, 'single_digit_digits.pdf')
dev.off()


#### Back up

# t-SNE of entire dataset
#Xg <- fftRtsne(PCs, 2,  df=1.0 ,rand_seed=4);
#p <- ggplot(data.frame(Xg,col=factor(y)), aes(x=X1,y=X2,color=col)) + geom_point(show.legend = FALSE)  +  theme(plot.title = element_text(size=8))

# t-SNE of entire dataset, showing only zeros
#Xg0.1 <- fftRtsne(PCs, 2,  df=0.1 ,rand_seed=4);
#idx <- 1:length(y)
#idx <- (y == 0)
#p <- ggplot(data.frame(Xg0.1[idx,],col=factor(y[idx])), aes(x=X1,y=X2,color=col)) + geom_point()  +  theme(plot.title = element_text(size=8))
#
#
#dfs <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
#
#idx <- 1:length(y)
#idx <- (y == 0 | y==3)
#
#idx <- 1:length(y)
##dfs <- c(0.01, 0.1,  0.25, 0.5, 0.75, 1.0)
#dfs <- c(-0.9,  -0.5, -0.25, 0,  0.5,  1.0 )
#
##dfs <- c(0.1)
#dfs <- c(-0.9,  -0.5, -0.25, 0,  0.5,  1.0 )
#dfs <- c(-0.999,-0.99,-0.95, -0.9)
#idx <- (y == 0)
#p_2 <- list()
#for(i in 1:length(dfs)){
#    #idx <- 1:length(y)
#    Xg_0 <- fftRtsne(PCs[idx,], 2,  df=dfs[i] ,rand_seed=4, max_iter=4000);
#    #idx <- (y == 0)
#   p_2[[i]] <-  ggplot(data.frame(Xg_0,col=factor(y[idx])), aes(x=X1,y=X2,color=col)) + geom_point(size=0.5)  +  theme(plot.title = element_text(size=8)) + ggtitle(sprintf("df=%f", dfs[i]))
#}
#do.call(grid.arrange,c(p_2,nrow=2))

