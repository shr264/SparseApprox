library(MASS)
source("projecthelper.r")
size <- 2^8
haarmat <- read.table("haarmat.txt",header = FALSE)
invhaarmat <- t(haarmat)
wvlt <- invhaarmat%*%data
write.table(matrix(haarmat%*%wvlt,nrow = size, ncol = size), file = "original.txt", row.names = FALSE, col.names = FALSE)
noisywvlt <- wvlt + rnorm(n=size^2,0,0.1)
write.table(matrix(haarmat%*%noisywvlt,nrow = size, ncol = size), file = paste(toString(0.1),"noisy.txt",sep=""), row.names = FALSE, col.names = FALSE)
recover <- admm(noisywvlt, invhaarmat, 0.001, rep(0,length=size^2), eps = 10^(-6), 1, data)
write.table(matrix(recover$betahat, nrow = size, ncol = size), file = paste(toString(0.1),"recovered.txt",sep=""), row.names = FALSE, col.names = FALSE)
