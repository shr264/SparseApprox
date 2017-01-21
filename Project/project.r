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

set.seed(12345)
p <- 200
N <- 1000
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

eta <- 0.1
initialt <- 1
subgradientdim2(y,X,eta,betaold,lambda,eps,initialt,betastar)

ista2a <- ista2(y,X,lambda,betaold,eps,betastar)

fista2a <- fista2(y,X,lambda,betaold,eps,betastar)

initialL <- 1
eta <- 0.5

istadim2(y,X,lambda,betaold, eps, initialL ,eta,betastar)

fistadim2(y,X,lambda,betaold, eps, initialL ,eta,betastar)

rho <- 1
admm2a <- admm2(y, X, lambda, betaold, eps = 10^(-3), rho, betastar)



cv5sg <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- subgradientconst2(y[-folds[[i]]],X[-folds[[i]],],betaold,lambda,eps,0.0001,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

cv5ista <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- ista2(y[-folds[[i]]],X[-folds[[i]],],lambda,betaold,eps,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

cv5fista <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- fista2(y[-folds[[i]]],X[-folds[[i]],],lambda,betaold,eps,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}


cv5admm <- function(y,X,betaold,eps,betastar,rho, nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            rho <- 1
            (fit <- admm(y[-folds[[i]]],X[-folds[[i]],],lambda, betaold, eps = 10^(-3), rho, betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

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

pdf('cvgc.pdf')
plot(ista2a$error,ylab=expression(group("||",beta^k - beta^(k-1),"||")[infinity]),xlab = 'Iteration',type = 'l', col = '1', ylim = c(0,0.3))
lines(fista2a$error, col = '2')
lines(subg2$error, col = '3')
lines(admm2a$error, col = '4')
legend('topright', c('ISTA', 'FISTA', 'SG','ADMM'), lty = c(1,1,1,1), col = c(1,2,3,4))
dev.off()

op <- microbenchmark(
    SG = subgradientconst2(y,X,betaold,lambdasg,eps,0.0001,betastar),
    ISTA = ista2(y,X,lambdaista,betaold,eps,betastar),
    FISTA = fista2(y,X,lambdafista,betaold,eps,betastar),
    ADMM = admm2(y, X, lambdaadmm, betaold, eps = 10^(-3), 1, betastar),
    times=10L)
print(op)
pdf('timing.pdf')
boxplot(op)
dev.off()

size <- 2^5
cameraman <- channel(readImage("cameraman.png"),"gray")
cameraman <- resize(cameraman,size,size)
data <- imageData(cameraman)
data <- matrix(data,nrow = size^2, ncol = 1)
haarmat <- GenW(size^2, family="DaubExPhase")
invhaarmat <- t(haarmat)
wvlt <- invhaarmat%*%data
min(data)
####
###savemat(matrix(invhaarmat%*%wvlt,nrow = size, ncol = size), filename = 'originalimage.jpeg', dev = 'jpeg')
writePNG(matrix(haarmat%*%wvlt,nrow = size, ncol = size), target = 'originalimage.png')
###addnoise
noisywvlt <- wvlt + rnorm(n=size^2,0,0.01)
min(noisywvlt)
###noisy image
###savemat(raster(matrix(invhaarmat%*%noisywvlt,nrow = size, ncol = size)),filename = 'noisy.jpeg', dev = 'jpeg')
writePNG(matrix(haarmat%*%noisywvlt,nrow = size, ncol = size), target = 'noisy.png')
### recovered image
recover <- admm(noisywvlt, invhaarmat, 0.01, rep(1,length=size^2), eps = 10^(-6), 1, data)
###savemat(matrix(recover$betahat, nrow = size, ncol = size),filename = 'recovered.jpeg', dev = 'jpeg')
writePNG(matrix(recover$betahat, nrow = size, ncol = size), target = 'recovered.png')

psnr(data,recover$betahat)

psnr(data,haarmat%*%noisywvlt)


library(clusterGeneration)
library(Matrix)
library(MASS)
p <- 10
spars <- 0.8
set.seed(12345)
Omega <- genPositiveDefMat(p)$Sigma
threshold <- quantile(as.vector(Omega), c(spars)) 
Omega[abs(Omega)<threshold] <- 0
Omega <- nearPD(Omega)$mat
Y <- mvrnorm(n = 2*p, mu = rep(0,p), Sigma = solve(Omega), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
 

xold <- diag(p)
rho <- 1
eps <- 10^(-3)
t <- 0.1

glassoista <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    xnew <- xold + t*(solve(xold)-S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    err <- max(abs(xnew-xold))
    itr <- 1
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        xnew <- xold + t*(solve(xold)-S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        err <- append(err,max(abs(xnew-xold)))
    }
return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassoista(Y,xold,0.1,lambda,eps=10^(-3),Omega)

    
eigen(Omega)



