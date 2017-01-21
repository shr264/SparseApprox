setwd("/Users/syedrahman/Documents/Fall2016/SparseApprox/Project")
library(MASS)
library(lassoshooting)
library(microbenchmark)
library("EBImage")
library("wavelets")
library("spatstat")
library("wavethresh")
library("squash")
library("raster")
library("png")
source("projecthelper.r")

p <- 200
for(N in c(500)){
set.seed(12345)
X <- mvrnorm(n=N,mu=rep(0,p),Sigma = diag(p))
indx <- sample(1:p,floor(p/10))
betastar <-  rep(0,p)
betastar[indx] <- runif(length(indx),5,50)
E <- rnorm(n=dim(X)[1],0,1)
y <- X%*%betastar + E

lambda <- 0.1
eps <- 10^(-3)
betaold <- rep(1,p)
t <- 0.01

initialt <- 1
#subgradientdim2(y,X,eta,betaold,lambda,eps,initialt,betastar)

#ista2a <- ista2(y,X,lambda,betaold,eps,betastar)

#fista2a <- fista2(y,X,lambda,betaold,eps,betastar)

initialL <- 1
eta <- 0.5

#istadim2(y,X,lambda,betaold, eps, initialL ,eta,betastar)

#fistadim2(y,X,lambda,betaold, eps, initialL ,eta,betastar)

rho <- 1
#admm2a <- admm2(y, X, lambda, betaold, eps = 10^(-3), rho, betastar)

nlambda <- seq(from = 0.0001, to = 1, length = 10)
lambdasg <- cv5sg(y,X,betaold,eps,betastar, nlambda)
lambdaista <- cv5ista(y,X,betaold,eps,betastar, nlambda)
lambdafista <- cv5fista(y,X,betaold,eps,betastar, nlambda)
lambdaadmm <- cv5admm(y,X,betaold,eps,betastar,1, nlambda)

subg2 <- subgradientconst2(y,X,betaold,lambdasg,eps,0.0001,betastar)
ista2a <- ista2(y,X,lambdaista,betaold,eps,betastar)
fista2a <- fista2(y,X,lambdafista,betaold,eps,betastar)
rho <- 1
admm2a <- admm2(y, X, lambdaadmm, betaold, eps = 10^(-3), rho, betastar)

pdf(paste(toString(N),'cvgc.pdf',sep=''))
plot(ista2a$error,ylab=expression(group("||",beta^k - beta^(k-1),"||")[infinity]),xlab = 'Iteration',type = 'l', col = '1', ylim = c(0,0.3))
lines(fista2a$error, col = '2')
lines(subg2$error, col = '3')
lines(admm2a$error, col = '4')
legend('topright', c('ISTA', 'FISTA', 'SG','ADMM'), lty = c(1,1,1,1), col = c(1,2,3,4))
dev.off()

pdf(paste(toString(N),'cvgc2.pdf',sep=''))
plot(ista2a$error2/sum(betastar^2),ylab='Relative Norm Error',xlab = 'Iteration',type = 'l', col = '1', ylim = c(0,0.0085))
lines(fista2a$error2/sum(betastar^2), col = '2')
lines(subg2$error2/sum(betastar^2), col = '3')
lines(admm2a$error2/sum(betastar^2), col = '4')
legend('topright', c('ISTA', 'FISTA', 'SG','ADMM'), lty = c(1,1,1,1), col = c(1,2,3,4))
dev.off()

op <- microbenchmark(
    SG = subgradientconst2(y,X,betaold,lambdasg,eps,0.0001,betastar),
    ISTA = ista2(y,X,lambdaista,betaold,eps,betastar),
    FISTA = fista2(y,X,lambdafista,betaold,eps,betastar),
    ADMM = admm2(y, X, lambdaadmm, betaold, eps = 10^(-3), 1, betastar),
    times=10L)
print(op)
pdf(paste(toString(N),'timing.pdf',sep=''))
boxplot(op)
dev.off()

difftable <- rbind(subg2$diff,ista2a$diff,fista2a$diff,admm2a$diff)/sum(betastar^2)
row.names(difftable) <- c( 'SG','ISTA', 'FISTA','ADMM')
write.table(difftable,file=paste(toString(N),'relativenormerror.txt',sep=''),row.names = TRUE, col.names = FALSE) 
}



