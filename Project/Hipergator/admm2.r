library(MASS)
source("projecthelper.r")
size <- 2^8
haarmat <- read.table("haarmat.txt",header = FALSE)
invhaarmat <- t(haarmat)
wvlt <- haarmat%*%data
write.table(matrix(invhaarmat%*%wvlt,nrow = size, ncol = size), file = "original.txt", row.names = FALSE, col.names = FALSE)
noisywvlt <- wvlt + rnorm(n=size^2,0,AA)
write.table(matrix(invhaarmat%*%noisywvlt,nrow = size, ncol = size), file = paste(toString(AA),"noisy.txt",sep=""), row.names = FALSE, col.names = FALSE)
recover <- admm(noisywvlt, haarmat, BB, rep(0,length=size^2), eps = 10^(-6), 1, data)
write.table(matrix(recover$betahat, nrow = size, ncol = size), file = paste(toString(AA),"recovered.txt",sep=""), row.names = FALSE, col.names = FALSE)